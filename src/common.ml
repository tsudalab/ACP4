(* Copyright (C) 2022, Francois Berenger
   Tsuda laboratory, The University of Tokyo,
   5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan. *)

open Printf

module A = BatArray
module BA = Bigarray
module BA1 = BA.Array1
module CLI = Minicli.CLI
module Fn = Filename
module Ht = BatHashtbl
module IntMap = BatMap.Int
module L = BatList
module LO = Line_oriented
module Log = Dolog.Log
module S = BatString

module SL = struct
  type t = bool * float (* (label, score) *)
  let create label score =
    (label, score)
  let get_score (_, s) = s
  let get_label (l, _) = l
  (* to do a decreasing sort of a score labels list *)
  let high_score_first_cmp (_, s1) (_, s2) =
    BatFloat.compare s2 s1
  let to_string (label, score) =
    sprintf "%d %f" label score
end

module ROC = Cpm.MakeROC.Make(SL)

let performance_metrics maybe_out_fn db_fn score_labels d_cut dx =
  ROC.rank_order_by_score_a score_labels;
  let auc = ROC.fast_auc_a score_labels in
  let ef_05 = ROC.fast_enrichment_factor 0.05 score_labels in
  let pm_05 = ROC.fast_power_metric_a 0.05 score_labels in
  let log_title =
    sprintf "AUC=%.3f EF5=%.2f PM5=%.2f DB=%s d_c=%g dx=%g"
      auc ef_05 pm_05 db_fn d_cut dx in
  let plot_title =
    sprintf "AUC=%.3f EF5=%.2f PM5=%.2f\\nDB=%s" auc ef_05 pm_05 db_fn in
  Log.info "%s" log_title;
  BatOption.may (fun out_fn ->
      Log.info "writing ROC to %s" out_fn;
      LO.with_out_file out_fn (fun out ->
          let roc = ROC.fast_roc_curve_a score_labels in
          A.iter (fun (fpr, tpr) ->
              fprintf out "%g %g\n" fpr tpr
            ) roc
        );
      let nb_actives = A.count_matching fst score_labels in
      let nb_decoys = (A.length score_labels) - nb_actives in
      Gnuplot.roc_curve plot_title out_fn nb_actives nb_decoys
    ) maybe_out_fn

let int_of_bool = function
  | true -> 1
  | false -> 0

type int32_array   = (int32, BA.int32_elt,   BA.c_layout) BA1.t
type float32_array = (float, BA.float32_elt, BA.c_layout) BA1.t

(* indexing ALDH1_25conf w/ One_band consumes about 15% of 32GM RAM
   when using OCaml native integers and floats (not 32b like now) *)
type sparse_float_fp = { indexes: int32_array;
                         values: float32_array }

let is_active name =
  S.starts_with name "active"

let output_code output count nb_dx (name, encoded) =
  let label_int = if is_active name then +1 else -1 in
  (* classification label expected by liblinear is
     +1 (active) or -1 (inactive) *)
  fprintf output "%+d" label_int;
  A.iteri (fun i_chan channel ->
      IntMap.iter (fun i_dx feat ->
          (* liblinear wants indexes to start at 1 --> 1 + ... *)
          let feat_idx = 1 + i_dx + (i_chan * nb_dx) in
          fprintf output " %d:%g" feat_idx feat
        ) channel
    ) encoded;
  (* terminate this molecule vector *)
  fprintf output "\n";
  incr count

(* like output_code, but to a newline-terminated string *)
let sprintf_code nb_dx (name, encoded): string =
  let buff = Buffer.create 1024 in
  let label_int = if is_active name then +1 else -1 in
  (* classification label expected by liblinear is
     +1 (active) or -1 (inactive) *)
  bprintf buff "%+d" label_int;
  A.iteri (fun i_chan channel ->
      IntMap.iter (fun i_dx feat ->
          (* liblinear wants indexes to start at 1 --> 1 + ... *)
          let feat_idx = 1 + i_dx + (i_chan * nb_dx) in
          bprintf buff " %d:%g" feat_idx feat
        ) channel
    ) encoded;
  (* terminate this molecule's vector *)
  Buffer.add_char buff '\n';
  Buffer.contents buff

let int32_array_of_list (l: int list): int32_array =
  let n = L.length l in
  let arr = BA1.create BA.Int32 BA.C_layout n in
  L.iteri (fun i x ->
      BA1.unsafe_set arr i (Int32.of_int x)
    ) l;
  arr

let float32_array_of_list (l: float list): float32_array =
  let n = L.length l in
  let arr = BA1.create BA.Float32 BA.C_layout n in
  L.iteri (fun i x ->
      BA1.unsafe_set arr i x
    ) l;
  arr

(* compute the average of an array of fps *)
let average_many max_dim fp_arr =
  let sum = A.create max_dim 0.0 in
  let num_fps = float (A.length fp_arr) in
  A.iter (fun fp ->
      let n = BA1.dim fp.indexes in
      assert(n = BA1.dim fp.values);
      for i = 0 to n - 1 do
        let j = Int32.to_int (BA1.unsafe_get fp.indexes i) in
        let x = BA1.unsafe_get fp.values i in
        sum.(j) <- sum.(j) +. x
      done
    ) fp_arr;
  (* average *)  
  let idxs = ref [] in
  let avg_vals = ref [] in
  for i = 0 to max_dim - 1 do
    (* prepare for sparse vector *)
    if sum.(i) > 0.0 then
      begin
        let avg = sum.(i) /. num_fps in
        idxs := i :: !idxs;
        avg_vals := avg :: !avg_vals
      end
  done;
  (* convert back to sparse_float_fp *)
  { indexes = int32_array_of_list   (L.rev !idxs) ;
    values  = float32_array_of_list (L.rev !avg_vals) }

(* transform the array of (float IntMap.t) into something more efficient
   to compute Tanimoto *)
let freeze nb_dx (name, encoded): (string * sparse_float_fp) =
  let idxs = ref [] in
  let vals = ref [] in
  A.iteri (fun i_chan channel ->
      IntMap.iter (fun i_dx feat ->
          (* liblinear wants indexes to start at 1 --> 1 + ... *)
          let feat_idx = 1 + i_dx + (i_chan * nb_dx) in
          idxs := feat_idx :: !idxs;
          vals := feat :: !vals
        ) channel
    ) encoded;
  (name, { indexes = int32_array_of_list   (L.rev !idxs) ;
           values  = float32_array_of_list (L.rev !vals) })

let tanimoto' fp1 fp2 =
  let icard = ref 0.0 in
  let ucard = ref 0.0 in
  let len1 = BA1.dim fp1.indexes in
  let len2 = BA1.dim fp2.indexes in
  let i = ref 0 in
  let j = ref 0 in
  while !i < len1 && !j < len2 do
    (* unsafe *)
    let k1 = BA1.unsafe_get fp1.indexes !i in
    let v1 = BA1.unsafe_get fp1.values  !i in
    let k2 = BA1.unsafe_get fp2.indexes !j in
    let v2 = BA1.unsafe_get fp2.values  !j in
    (* process keys in increasing order *)
    if k1 < k2 then
      (ucard := !ucard +. v1;
       incr i)
    else if k2 < k1 then
      (ucard := !ucard +. v2;
       incr j)
    else (* k1 = k2 *)
    if v1 <= v2 then
      (icard := !icard +. v1;
       ucard := !ucard +. v2;
       incr i;
       incr j)
    else
      (icard := !icard +. v2;
       ucard := !ucard +. v1;
       incr i;
       incr j)
  done;
  while !i < len1 do (* finish fp1; unsafe *)
    ucard := !ucard +. (BA1.unsafe_get fp1.values !i);
    incr i
  done;
  while !j < len2 do (* finish fp2; unsafe *)
    ucard := !ucard +. (BA1.unsafe_get fp2.values !j);
    incr j
  done;
  if !icard = 0.0 then
    0.0 (* NaN protection: ucard=0 --> icard=0 *)
  else
    !icard /. !ucard

let tanimoto (_l1, fp1) (_l2, fp2) =
  tanimoto' fp1 fp2

let tani_dist' x y =
  1.0 -. (tanimoto' x y)

let tani_dist x y =
  1.0 -. (tanimoto x y)

let parse_one verbose cutoff dx nb_dx input: string * sparse_float_fp =
  let cand_ph4_mol = Ph4.read_one_ph4_encoded_molecule input in
  let encoded = Ph4.encode verbose cutoff dx cand_ph4_mol in
  freeze nb_dx encoded

let parse_all verbose cutoff dx nb_dx input_fn =
  LO.with_in_file input_fn (fun input ->
      let res, exn =
        L.unfold_exn (fun () ->
            parse_one verbose cutoff dx nb_dx input
          ) in
      assert(exn = End_of_file);
      res
    )

let read_one_bin input: string * sparse_float_fp =
  Marshal.from_channel input

let write_one_bin output (x: string * sparse_float_fp): unit =
  Marshal.(to_channel output x [No_sharing])

module Ligand_defaults = struct
  let radial_cutoff = 5.0
  let dx = 0.5
  let bst_chunk_size = 100_000
end

(* protein Binding-Sites *)
module BS_defaults = struct
  let radial_cutoff = 40.0
  let dx = 0.9
end

let stddev (a: float array): float =
  let sqr x =
    x *. x in
  let n, sx, sx2 =
    A.fold_left
      (fun (n, sx, sx2) x -> (succ n, sx +. x, sx2 +. (sqr x)))
      (0, 0., 0.) a
  in
  sqrt ((sx2 -. (sqr sx) /. (float n)) /. (float n))
