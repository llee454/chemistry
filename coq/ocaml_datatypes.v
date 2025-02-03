Require Coq.extraction.Extraction.
Require Import Reals.
Require Import List.
Require Import Strings.String.

Open Scope list_scope.

Extract Inductive prod => "(*)"  [ "(,)" ].

Extract Constant fst => "fst".

Extract Constant snd => "snd".

Extraction Inline fst snd.

Extract Inductive list => "list" ["[]" "(fun (x, xs) -> List.cons x xs)"]
  "(fun branch0 branch1 xs ->
    match xs with
    | [] -> branch0 ()
    | x :: xs -> branch1 x xs)".

Extract Constant app => "List.append".

Extraction Inline app.

Extract Inductive sumbool => "bool" [ "true" "false" ].

Extract Inductive bool => "bool" ["true" "false"]
  "(fun branch0 branch1 b ->
     match b with
     | true -> branch0 ()
     | false -> branch1 ())".

Extract Constant negb => "Bool.not".

Extract Constant String.eqb => "String.equal".

Extract Inductive string => "string" ["" "UNDEFINED_TERM"].

Extract Constant existsb => "(fun f xs -> Core_kernel.List.exists xs ~f)".

Extract Constant filter => "(fun f xs -> Core_kernel.List.filter xs ~f)".

Module Dummy_reals <: RbaseSymbolsSig.

  Declare Scope R_scope.

  Delimit Scope R_scope with R.

  Inductive R_ : Set := rdummy.

  Extract Inductive R_ => "float" [ "0.0" ].

  Definition of_number (_ : Number.int) := rdummy.

  Definition to_number (_ : R_) := Number.IntDecimal (Decimal.Pos (Decimal.D0 Decimal.Nil)).

  Number Notation R_ of_number to_number : R_scope.

  Definition R := R_.

  Extract Inlined Constant R => "float".

  Axiom Rabst : ConstructiveCauchyReals.CReal -> R.

  (* Definition Rrepr (_ : R) := ConstructiveCauchyReals.mkCReal (fun (_ : Z) => QArith_base.inject_Z Z0) Z0. *)

  Axiom Rrepr :  R -> ConstructiveCauchyReals.CReal.

  Axiom Rquot1 : forall x y : R, ConstructiveCauchyReals.CRealEq (Rrepr x) (Rrepr y) -> x = y.

  Axiom Rquot2 : forall x : ConstructiveCauchyReals.CReal, ConstructiveCauchyReals.CRealEq (Rrepr (Rabst x)) x.

  Definition R0 : R := rdummy.

  Extract Inlined Constant R0 => "0.0".

  Definition R1 : R := rdummy.

  Extract Inlined Constant R1 => "1.0".

  Definition Rplus (_ _: R) : R := rdummy.

  Extract Inlined Constant Rplus => "( +. )".

  Definition Rmult (_ _ : R) : R := rdummy.

  Extract Inlined Constant Rmult => "( *. )".

  Definition Ropp (_ : R) : R := rdummy.

  Extract Inlined Constant Ropp => "-".

  Definition Rlt (_ _ : R) : Prop := False.

  Extract Inlined Constant Rlt => "( <. )".

  Axiom R0_def : R0 = Rabst (ConstructiveCauchyReals.inject_Q {| QArith_base.Qnum := 0; QArith_base.Qden := 1 |}).

  Axiom R1_def : R1 = Rabst (ConstructiveCauchyReals.inject_Q {| QArith_base.Qnum := 1; QArith_base.Qden := 1 |}).

  Axiom Rplus_def : forall x y : R, Rplus x y = Rabst (ConstructiveCauchyReals.CReal_plus (Rrepr x) (Rrepr y)).

  Axiom Rmult_def : forall x y : R, Rmult x y = Rabst (ConstructiveCauchyRealsMult.CReal_mult (Rrepr x) (Rrepr y)).

  Axiom Ropp_def : forall x : R, Ropp x = Rabst (ConstructiveCauchyReals.CReal_opp (Rrepr x)).

  Axiom Rlt_def : forall x y : R, Rlt x y = ConstructiveCauchyReals.CRealLtProp (Rrepr x) (Rrepr y).

  Infix "+" := Rplus : R_scope.

  Infix "*" := Rmult : R_scope.

  Notation "- x" := (Ropp x) : R_scope.

  Infix "<" := Rlt : R_scope.

End Dummy_reals.