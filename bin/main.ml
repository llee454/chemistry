open! Core
open! Chemistry

let () = Bohr.shielding_constant ~z:2. ~n:1. (electron_volts_to_joules 24.58) |> printf "%0.4f\n"
