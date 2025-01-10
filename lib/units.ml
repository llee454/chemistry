(**
  This module illustrates how we can use GADTs to enforce proper
  units.
*)
open! Core

module Measurement = struct
  type kg = Kg

  type g = G

  (** Meters per second squared *)
  type mps2 = Mps2

  (** Newtons *)
  type n = N

  type _ t =
    | Kg : float -> kg t
    | G : float -> g t
    | Mps2 : float -> mps2 t
    | N : float -> n t
end

module Use_units = struct
  open! Measurement

  let get_force (Kg mass) (Mps2 acceleration) : n t = N (mass *. acceleration)
end
