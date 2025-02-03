(**
  Defines a set of functions for performing chemistry calculations.
*)
Require Import Coq.Strings.String.
Require Import List.
Require Import micromega.Lra.
Import ListNotations.
Require Coq.extraction.Extraction.

Require Import ocaml_datatypes.

Open Scope list_scope.

Module Chemistry.

Include Dummy_reals.

Open Scope R_scope.

Inductive unit : Set := Meter | Angstrom.

Definition measurement_unit : Set := list unit * list unit.

Module Measurement.

  Inductive t : measurement_unit -> Set :=
  | Cast : forall (_ : R) (units : measurement_unit), t units.

  Definition value {units : measurement_unit} (x : t units) :=
    match x with
    | Cast value _ => value
    end.
  
End Measurement.

Definition angstrom_to_meters_const : R := 10_000_000_000%R.

Extract Constant angstrom_to_meters_const => "1e10".

Definition angstrom_to_meters (x : Measurement.t ([Angstrom], nil)) : R :=
  (Measurement.value x) * angstrom_to_meters_const.

End Chemistry.