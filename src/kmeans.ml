(* Copyright (C) 2023, Francois Berenger
   Tsuda laboratory, The University of Tokyo,
   5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.

   k-means clustering and silhouette metric.

   Silhouettes: A graphical aid to the interpretation and validation
   of cluster analysis
   https://doi.org/10.1016/0377-0427(87)90125-7 *)

open Printf

module A = BatArray
module CLI = Minicli.CLI
module Ht = BatHashtbl
module IntSet = BatSet.Int
module L = BatList
module LO = Line_oriented
module Log = Dolog.Log

type 'a cluster = { members: IntSet.t; (* indexes of all cluster members *)
                    center: 'a } (* the cluster center (average in k-means) *)

let get_center c =
  c.center

let get_members c =
  c.members

(* randomly choose [k] elements from [arr] *)
let rand_choose rng k arr =
  let shuffled = A.copy arr in
  A.shuffle ~state:rng shuffled;
  A.left shuffled k

let rand_select rng arr =
  let n = A.length arr in
  let i = Random.State.int rng n in
  arr.(i)

(* return cluster id for x (nearest center's index) *)
let assign_to_cluster rng dist centers x =
  let distances = A.map (dist x) centers in
  let mini = A.min distances in
  let res = ref [] in
  A.iteri (fun i d ->
      if d = mini then
        res := i :: !res
    ) distances;
  let candidates = A.of_list !res in
  if A.length candidates > 1 then
    (Log.warn "Kmeans.assign_to_cluster: equidistant";
     rand_select rng candidates)
  else
    candidates.(0)

(* from an assignment (a mapping from elt. index to cluster_id), create
   a list of clusters *)
let extract_clusters assignment =
  let cid2cluster = Ht.create 11 in
  A.iteri (fun elt_i cid ->
      let prev = Ht.find_default cid2cluster cid IntSet.empty in
      let curr = IntSet.add elt_i prev in
      Ht.replace cid2cluster cid curr
    ) assignment;
  L.map snd (Ht.bindings cid2cluster)

let compute_center max_dim all_elements cluster_members =
  let members = IntSet.to_array cluster_members in
  let element_members = A.map (fun i -> all_elements.(i)) members in
  Common.average_many max_dim element_members

let cluster_variance dist all_elements cluster =
  let n = IntSet.cardinal cluster.members in
  let sum_d2 =
    IntSet.fold (fun i acc ->
        let d = dist cluster.center all_elements.(i) in
        acc +. (d *. d)
      ) cluster.members 0.0 in
  sum_d2 /. (float n)

let log_clusters clusts vars =
  let i = ref 0 in
  L.iter2 (fun c v ->
      Log.info "c%d: n=%d var=%f" !i (IntSet.cardinal c.members) v;
      incr i
    ) clusts vars

let max_iter = 100
let epsilon = 0.001

(* automate search for optimal k *)
(* implement k-range: start:step:stop *)

(* FBR: examine visualy cluster quality for some targets *)
(* FBR: after optimal k has been determined on a dataset,
 *      compute the average hit-rate per cluster, for clusters
 *      with at-least one active; might relate to the difficulty
 *      of a dataset *)

(* the k-means implementation *)
let cluster max_dim rng dist k elements =
  let () =
    let n = A.length elements in
    if k > n then
      (Log.fatal "Kmeans.cluster: k=%d > n=%d" k n;
       exit 1) in
  (* random initial cluster centers *)
  let init_centers = rand_choose rng k elements in
  (* compute initial assignment *)
  let clusters =
    let assignment =
      A.map (assign_to_cluster rng dist init_centers) elements in
    extract_clusters assignment in
  let cluster_w_centers =
    L.map2 (fun members center -> { members; center })
      clusters (A.to_list init_centers) in
  let vars = L.map (cluster_variance dist elements) cluster_w_centers in
  Log.info "avg_var: %f" (L.favg vars);
  let rec loop iter clusts vars =
    if iter >= max_iter then
      (Log.info "Kmeans.cluster: max iter";
       (clusts, vars))
    else
      let centers' = A.of_list (L.map get_center clusts) in
      (* compute assignment *)
      let clusters' =
        let assignment =
          A.map (assign_to_cluster rng dist centers') elements in
        extract_clusters assignment in
      (* update centers *)
      let clusts' =
        let new_centers = L.map (compute_center max_dim elements) clusters' in
        L.map2 (fun members center -> { members; center })
          clusters' new_centers in
      let prev_var = L.favg vars in
      let vars' = L.map (cluster_variance dist elements) clusts' in
      let curr_var = L.favg vars' in
      Log.info "avg_var: %f" curr_var;
      let delta = prev_var -. curr_var in
      (* Log.info "Kmeans.cluster: delta = %f" delta; *)
      if abs_float delta < epsilon then
        (clusts', vars')
      else
        loop (succ iter) clusts' vars' in
  (* repeat until max_iter or convergence *)
  loop 0 cluster_w_centers vars

(* average distance of elt_i to all _other_ cluster elements *)
let avg_dist dist cluster i =
  let res = ref [] in
  IntSet.iter (fun j ->
      if i <> j then
        let d = dist i j in
        res := d :: !res
    ) cluster.members;
  L.favg !res

(* the silhouette metric *)
let silhouette dist clusters i =
  let my_cluster', other_clusters =
    L.partition (fun c -> IntSet.mem i c.members) clusters in
  assert(L.length my_cluster' = 1);
  let my_cluster = L.hd my_cluster' in
  if IntSet.cardinal my_cluster.members = 1 then
    0.0 (* from the definition *)
  else
    (* a: average dist to its own cluster members *)
    let a = avg_dist dist my_cluster i in
    (* b: min(avg dist to other clusters) *)
    let avg_dists =
      L.map (fun c ->
          avg_dist dist c i
        ) other_clusters in
    let b = L.min avg_dists in
    (b -. a) /. (max a b)

let avg_silhouette dist clusters =
  let n = L.sum (L.map (fun c -> IntSet.cardinal c.members) clusters) in
  let silhouettes = A.init n (silhouette dist clusters) in
  A.favg silhouettes

let find_max_dim mols =
  let high_indexes = A.map Common.high_index mols in
  1 + (Int32.to_int (A.max high_indexes))

(* clusters textual output: '^mol_name\tcID$' *)
let write_clusters_out fn clusters all_names =
  LO.with_out_file fn (fun output ->
      L.iteri (fun i c ->
          IntSet.iter (fun j ->
              fprintf output "%s\tc%d\n" all_names.(j) i
            ) c.members
        ) clusters
    )

(* useful only for silhouette metric calculation; since *)
(* during clustering the distance of interest is to *)
(* the cluster center only (which is not a point from *)
(* the dataset) *)
let tani_cached all_mols mtx i j =
  let d = mtx.(i).(j) in
  if d > -1.0 then
    d (* already in cache *)
  else
    (* init cache *)
    let d' = Common.tani_dist' all_mols.(i) all_mols.(j) in
    mtx.(i).(j) <- d';
    mtx.(j).(i) <- d';
    d'

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n\
              %s\n  \
              -i <filename.ph4>: input file\n  \
              -o <filename.txt>: output file\n  \
              [-np <int>]: maximum number of CPU cores (default=1)\n  \
              [-k int]: specify the number of clusters\n  \
              [--ks 'start:step:stop']: scan range for best k\n  \
              [-v]: verbose/debug mode\n"
       Sys.argv.(0);
     exit 1);
  let input_fn = CLI.get_string ["-i"] args in
  let output_fn = CLI.get_string ["-o"] args in
  let _nprocs = CLI.get_int_def ["-np"] args 1 in
  let k = CLI.get_int ["-k"] args in
  let _k_range_str = CLI.get_string_opt ["--ks"] args in
  let verbose = CLI.get_set_bool ["-v"] args in
  let binding_site_mode = CLI.get_set_bool ["--BS"] args in
  let cutoff =
    CLI.get_float_def ["-c"] args
      (if binding_site_mode
       then Common.BS_defaults.radial_cutoff
       else Common.Ligand_defaults.radial_cutoff) in
  let dx =
    CLI.get_float_def ["-dx"] args
      (if binding_site_mode
       then Common.BS_defaults.dx
       else Common.Ligand_defaults.dx) in
  let rng = match CLI.get_int_opt ["-s"] args with
    | None -> Random.State.make_self_init ()
    | Some seed -> Random.State.make [|seed|] in
  CLI.finalize (); (* ------------------------------------------------------ *)
  let nb_dx = 1 + BatFloat.round_to_int (cutoff /. dx) in
  let max_dim = 1 + (Ph4.nb_channels * nb_dx) in
  Log.info "reading molecules...";
  let names, all_mols =
    A.split (A.of_list (Common.parse_all verbose cutoff dx nb_dx input_fn)) in
  let nb_mols = A.length all_mols in
  let dist_cache = A.make_matrix nb_mols nb_mols (-1.0) in
  let dist = tani_cached all_mols dist_cache in
  let act_max_dim = find_max_dim all_mols in
  Log.info "max_dim: %d act_max_dim: %d" max_dim act_max_dim;
  let clusts, vars = cluster act_max_dim rng Common.tani_dist' k all_mols in
  log_clusters clusts vars;
  write_clusters_out output_fn clusts names;
  let sil = avg_silhouette dist clusts in
  Log.info "avg_sil: %f" sil

let () = main ()
