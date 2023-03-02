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

(* randomly choose [k] elements from [arr] *)
let rand_choose rng k arr =
  let shuffled = A.copy arr in
  A.shuffle ~state:rng shuffled;
  A.left shuffled k

(* return cluster id for x (first nearest center's index) *)
let assign_to_cluster dist centers x =
  let mini = ref infinity in
  let j = ref (-1) in
  let distances = A.map (dist x) centers in
  A.iteri (fun i d ->
      if d < !mini then
        (j := i;
         mini := d)
    ) distances;
  !j

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

let euclid xs ys =
  let sum_diff2 =
    L.fold_left2 (fun acc x y ->
        let d = x -. y in
        acc +. (d *. d)
      ) 0.0 xs ys in
  sqrt sum_diff2

(* FBR: log out the clusters: *)
(*      cid, members, variance, silhouette *)

(* the k-means implementation *)
let cluster max_dim rng dist k elements =
  (* random initial cluster centers *)
  let init_centers = rand_choose rng k elements in
  (* compute initial assignment *)
  let clusters =
    let assignment = A.map (assign_to_cluster dist init_centers) elements in
    extract_clusters assignment in
  (* compute centers *)
  let centers = L.map (compute_center max_dim elements) clusters in
  let cluster_w_centers =
    L.map2 (fun members center -> { members; center }) clusters centers in
  let variances = L.map (cluster_variance dist elements) cluster_w_centers in
  let rec loop iter clusts vars =
    if iter >= max_iter then
      (Log.info "Kmeans.cluster: max iter";
       (clusts, vars))
    else
      let centers' = A.of_list (L.map (fun clust -> clust.center) clusts) in
      (* compute assignment *)
      let clusters' =
        let assignment = A.map (assign_to_cluster dist centers') elements in
        extract_clusters assignment in
      (* update centers *)
      let clusts' =
        let new_centers = L.map (compute_center max_dim elements) clusters' in
        L.map2 (fun members center -> { members; center })
          clusters' new_centers in
      let vars' = L.map (cluster_variance dist elements) clusts' in
      let delta = euclid vars vars' in
      Log.info "Kmeans.cluster: delta = %f" delta;
      if delta < epsilon then
        (clusts', vars')
      else
        loop (succ iter) clusts' vars' in
  (* repeat until max_iter or convergence *)
  loop 0 cluster_w_centers variances

let find_max_dim mols =
  let high_indexes = A.map Common.high_index mols in
  1 + (Int32.to_int (A.max high_indexes))

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n\
              %s\n  \
              -i <filename.ph4>: input file\n  \
              [-np <int>]: maximum number of CPU cores (default=1)\n  \
              [-k int]: specify the number of clusters\n  \
              [-v]: verbose/debug mode\n"
       Sys.argv.(0);
     exit 1);
  let input_fn = CLI.get_string ["-i"] args in
  let _nprocs = CLI.get_int_def ["-np"] args 1 in
  let k = CLI.get_int ["-k"] args in
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
  let _names, all_mols =
    A.split (A.of_list (Common.parse_all verbose cutoff dx nb_dx input_fn)) in
  let act_max_dim = find_max_dim all_mols in
  Log.info "max_dim: %d act_max_dim: %d" max_dim act_max_dim;
  let _nb_mols = A.length all_mols in
  let clusts, vars = cluster act_max_dim rng Common.tani_dist' k all_mols in
  log_clusters clusts vars

let () = main ()
