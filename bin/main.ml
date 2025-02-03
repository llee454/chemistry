open! Core
open! Chemistry
open! Float

(* let () = Bohr.shielding_constant ~z:2. ~n:1. (electron_volts_to_joules 24.58) |> printf "%0.4f\n" *)
(*
let () = 
  let wavelength = 5_890.0 in
  let energy = light_wavelength_to_frequncy wavelength |> photon_energy in
  printf !"photon energy: %fE-29 j\n" (energy * 1e29);
  let ionization_energy = electron_volts_to_joules 5.138 in
  printf !"ionization energy: %fE-19 j\n" (ionization_energy * 1e19);
  let shielding = Bohr.shielding_constant ~z:11.0 ~n:3.0 ionization_energy in
  let acc = ref None in
  for n = 1 to 10 do
    let excited_ionization_energy = Bohr.electron_energy ~z:(11.0 - shielding) (float n) in
    let energy_diff = abs (energy - (ionization_energy - excited_ionization_energy)) in
    printf !"n = %d excited_ionization energy = %fE-20J energy_diff = %fE-20J\n" n (excited_ionization_energy * 1e20) (energy_diff * 1e20);
    acc := match !acc with
    | None -> Some (n, energy_diff)
    | Some (prev_n, prev_energy_diff) as prev ->
      if energy_diff < prev_energy_diff
      then Some (n, energy_diff)
      else prev
  done;
  printf !"estimated excited orbit's principal quantum number: %d\n" (fst (Option.value_exn !acc))
*)
let () =
  printf !"Test: %f\n" Chemistry.(angstrom_to_meters (Measurement.Cast (1.0, ([Angstrom], []))));

