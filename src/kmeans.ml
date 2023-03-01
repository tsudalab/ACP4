(* Copyright (C) 2023, Francois Berenger
   Tsuda laboratory, The University of Tokyo,
   5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.

   k-means clustering and silhouette metric.

   Silhouettes: A graphical aid to the interpretation and validation of cluster analysis
   https://doi.org/10.1016/0377-0427(87)90125-7 *)

open Printf

module A = BatArray
module CLI = Minicli.CLI
module Ht = Hashtbl
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

let extract_clusters assignment =
  let cluster_ids = IntSet.of_array assignment in
  IntSet.fold (fun id acc ->
      let cluster =
        A.fold_lefti (fun acc i id' ->
            if id = id' then
              IntSet.add i acc
            else
              acc
          ) IntSet.empty assignment in
      cluster :: acc
    ) cluster_ids []

let compute_center add div cluster_members all_elements =
  let members = IntSet.to_array cluster_members in
  let n = A.length members in
  let element_members = A.map (fun i -> all_elements.(i)) members in
  let sum = A.reduce add element_members in
  div sum n

let max_iter = 100
let epsilon = 0.001

let cluster_variance dist cluster all_elements =
  let n = IntSet.cardinal cluster.members in
  let sum_d2 =
    IntSet.fold (fun i acc ->
        let d = dist cluster.center all_elements.(i) in
        acc +. (d *. d)
      ) cluster.members 0.0 in
  sum_d2 /. (float n)

let cluster rng add div dist k elements =
  let init_centers = rand_choose rng k elements in
  let assignment = A.map (assign_to_cluster dist init_centers) elements in
  let clusters = extract_clusters assignment in
  let _centers =
    L.map (fun cluster -> compute_center add div cluster elements) clusters in
  (* recompute centers *)
  (* recompute assignment *)
  (* as long as assignments are changing, refresh them *)
  failwith "not implemented yet"

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
  let nprocs = CLI.get_int_def ["-np"] args 1 in
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
  CLI.finalize (); (* ------------------------------------------------------ *)
  let nb_dx = 1 + BatFloat.round_to_int (cutoff /. dx) in
  Log.info "reading molecules...";
  let _names, all_mols =
    A.split (A.of_list (Common.parse_all verbose cutoff dx nb_dx input_fn)) in
  let nb_mols = A.length all_mols in
  let _all_indexes = A.init nb_mols (fun i -> i) in
  (* compute Gram matrix in // *)
  let matrix = A.make_matrix nb_mols nb_mols 0.0 in
  Log.info "Gram matrix initialization...";
  Molenc.Gram.initialize_matrix Common.tani_dist' nprocs 1 all_mols matrix;
  Molenc.Gram.print_corners matrix;
  Log.info "Adding nodes and edges to graph...";
  failwith "not implemented yet"

let () = main ()
