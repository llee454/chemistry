open! Core
open Float
 
include Coq_chemistry
 
(** The number of meters spanned by an angstrom *)
let angstrom = 1e-10

let angstroms_to_meters = ( * ) angstrom

let meters_to_angstroms x = x / angstrom

(** Planck's constant measured in Joule second *)
let plancks_constant = 0.66252e-33

let planck_bar = plancks_constant / (2. * pi)

(** The speed of light measured in meters per second. *)
let light_velocity = 299_792_458.

(** The avogrado's number of atoms. *)
let mole = 6.0229e23

(** The mass of a proton in kilograms. *)
let proton_mass = 1.67239e-27

(** The mass of a neutron in kilograms. *)
let neutron_mass = 1.67470e-27

(** The mass of a dalton in kilograms. *)
let dalton_to_kg = 1.66033e-27

module Constant = struct
  let h = plancks_constant

  let c = light_velocity

  let n = mole

  let d = dalton_to_kg
end

module Symbols = struct
  (** Note: this must be represented by a multibyte Unicode string. *)
  let angstrom_char = "\u{212B}"
end

let stoney_to_coulomb = 1.054822e-5

let columnb_to_stoney = 1. / stoney_to_coulomb

(** The mass of an electron in kg. *)
let electron_mass_kg = 0.91083e-30

(** The charge of an electron in coulombs. *)
let electron_charge_c = 0.160206e-18

let electron_charge_s = electron_charge_c * columnb_to_stoney

(** The number of joules represented by one electron volt *)
let electron_volt = 1.60206e-19

let electron_volts_to_joules = ( * ) electron_volt

let joules_to_electron_volts x = x / electron_volt

(**
  Calculates the potential energy (measured in newtons) between two
  charges [q0] and [q1] that are separated by [r] meters where [q0]
  and [q1] are measured in stoneys.
*)
let electrostatic_potential_energy_s ~q0 ~q1 ~r = q0 * q1 / r

(**
  Accepts the frequency of light (measured in cycles per second)
  and returns its wavelength measured in meters.
*)
let light_frequency_to_wavelength freq = light_velocity / freq

(**
  Accepts the wavelength (measured in meters) of light and returns
  its frequency (measured in cycles per second).
*)
let light_wavelength_to_frequncy wavelength = light_velocity / wavelength

(**
  Accepts the frequency of a photon (measured in cycles per second)
  and returns its energy measured in Joules.
*)
let photon_energy freq = plancks_constant * freq

let%expect_test "photon_energy" =
  1_216.0
  |> angstroms_to_meters
  |> light_wavelength_to_frequncy
  |> photon_energy
  |> joules_to_electron_volts
  |> printf "%0.4f eV";
  [%expect {| 10.1955 eV |}]

(**
  Accepts the velocity of an electron and returns its De Broglie
  wavelength.
*)
let de_broglie_wavelength velocity = plancks_constant / (electron_mass_kg * velocity)

let%expect_test "energy_of_photon" =
  let v = 650e-9 in
  light_wavelength_to_frequncy v |> photon_energy |> ( * ) 1e19 |> printf "%0.2fx10^19 J";
  [%expect {| 3.06x10^19 J |}]

let%expect_test "orbit_transition" =
  1_216. * angstrom
  |> light_wavelength_to_frequncy
  |> photon_energy
  |> ( * ) 1e18
  |> printf "%0.2fx10^18 J";
  [%expect {| 1.63x10^18 J |}]

(**
  Equations from the Bohr model.

  Note: the Bohr model is only approximately correct when applied
  to atoms that have a single electron.
*)
module Bohr = struct
  (**
    Accepts the principal quantum number (orbit) of an electron
    in a single electron.
  *)
  let angular_momentum n = n * planck_bar

  (**
    Accepts one argument, [nucleus_mass] the mass of the nucleus
    in kilograms, and returns the "reduced electron mass."
    
    Whenever, you use the electron mass in the Bohr Sommerfeld model,
    you can produce more accurate results by using the reduced
    electron mass, which is computed by accounting for the fact
    that both the electron and the proton orbit around a shared
    point between them (in the model).
  *)
  let reduced_electron_mass nucleus_mass = 1.0 / ((1.0 / electron_mass_kg) + (1.0 / nucleus_mass))

  (** The Bohr radius (a0) in meters *)
  let bohr_radius_constant ?(electron_mass = electron_mass_kg) () =
    square planck_bar / (electron_mass * square electron_charge_s)

  let%expect_test "bohr_radius" =
    printf "%0.4f %s" (bohr_radius_constant () / angstrom) Symbols.angstrom_char;
    [%expect {| 0.5292 Å |}]

  (**
    Accepts two arguments: z, the atomic number of an atom; and n,
    the principal quantum number (orbit) of an electron; and returns
    the radius of the orbit in meters.
  *)
  let bohr_radius ?electron_mass ~z n = square n / z * bohr_radius_constant ?electron_mass ()

  let%expect_test "bohr_radius" =
    bohr_radius ~z:(3.0 - 1.74) 2.0 |> meters_to_angstroms |> printf "%0.4f \u{212B}";
    [%expect {| 1.6799 Å |}]

  (**
    Accepts two arguments: [z], the atomic number of an atom;
    [n], the principal quantum number (orbit) of an electron;
    and returns the velocity of the electron in the given orbit
    measured in meters per second.
  *)
  let electron_velocity ~z n = z * square electron_charge_s / (n * planck_bar)

  let%expect_test "electron_velocity" =
    electron_velocity ~z:1.0 1.0 |> ( * ) 1e-6 |> printf "%0.4fx10^6";
    [%expect {| 2.1877x10^6 |}]

  (**
    Accepts three arguments: [z], the atomic number of an atom; [n],
    the principal quantum number (orbit) of an electron; and [l],
    the electron's orbit angular momentum quantum number; and
    returns the closest approach that the electron gets to the
    nucleus measured in meters.

    Note: Bohr, modeled the electron orbits using circular
    orbits. Sommerfield, extended this model to include eliptical
    orbits. This function uses a variant of the Bohr-Sommerfield
    model.

    Note: [z] can be written as [z - s] to reflect the idea of
    "shielding." In the Bohr-Sommerfield model, an electron can be
    partially shielded from the nucleus by other electrons that are
    in orbits closer to the nucleus than it. We subtract 1 from z
    for every electron that fully shields an outer electron. Certain
    outer orbits however can "penetrate" within inner orbits. When
    this happens the shielding constant s will be less than one.
  *)
  (* let closest_radius ~z ~n l =
       (square n) * bohr_radius_constant* (1.5 - (((square l) + l)/(2. * square n)))/z

     let%expect_test "closest_radius" =
       (closest_radius ~z:1. ~n:3. 2.) / angstrom |> printf !"%0.4f A\n";
       bohr_radius_constant / angstrom |> printf !"%0.4f";
       [%expect {||}] *)

  (** The energy of an electron in the 1s orbital of a hydrogen atom. *)
  let electron_energy_constant = electron_mass_kg * int_pow electron_charge_s 4 / (2. * square planck_bar)

  let%expect_test "electron_energy_constant" =
    printf "%0.4f eV" @@ joules_to_electron_volts electron_energy_constant;
    [%expect {| 13.6047 eV |}]

  (**
    Accepts two arguments: z, the atomic number of an atom; and n,
    the principal quantum number (orbit) of an electron; and returns
    the total energy of the given electron (measured in joules).

    Note: this equals the ionization energy needed to remove the
    electron.
  *)
  let electron_energy ~z n = -electron_energy_constant * square (z / n)

  (**
    Accepts two arguments: z, the atomic number of an atom; and n,
    the principal quantum number (orbit) of an electron; and returns
    the velocity (in meters per second) of the electron according
    to the Bohr model.
  *)
  let electron_velocity ~z n = z // n * (square electron_charge_s / planck_bar)

  let%expect_test "electron_velocity" =
    electron_velocity ~z:1 1 |> ( * ) 1e-6 |> printf "%0.4fx10^-6 m/s";
    [%expect {| 2.1877x10^-6 m/s |}]

  (**
    Accepts three arguments: [z], the atomic number of an atom;
    [n], the principal quantum number (orbit) of an electron; and
    [energy] the energy, measured in joules, needed to remove an
    electron from the atom (the ionization energy); and returns
    the shielding constant for the electron.
  *)
  let shielding_constant ~z ~n energy = z - (n * sqrt (energy / electron_energy_constant))

  let%expect_test "shielding_constant" =
    electron_volts_to_joules 24.58 |> shielding_constant ~z:2. ~n:1. |> printf "%0.4f\n";
    electron_volts_to_joules 2_085. |> shielding_constant ~z:13. ~n:1. |> printf "%0.4f";
    [%expect {|
      0.6559
      0.6203
      |}]
end

module Electronegativity = struct
  let tbl =
    [
      "H", 2.1;
      "Li", 1.0;
      "Be", 1.5;
      "B", 2.0;
      "C", 2.5;
      "N", 3.0;
      "O", 3.5;
      "F", 4.0;
      "Na", 0.9;
      "Mg", 1.2;
      "Al", 1.5;
      "Si", 1.8;
      "P", 2.1;
      "S", 2.5;
      "Cl", 3.0;
      "K", 0.8;
      "Ca", 1.0;
      "Sc", 1.3;
      "Ti", 1.5;
      "V", 1.6;
      "Cr", 1.6;
      "Mn", 1.5;
      "Fe", 1.8;
      "Co", 1.9;
      "Ni", 1.9;
      "Cu", 1.9;
      "Zn", 1.6;
      "Ga", 1.6;
      "Ge", 1.8;
      "As", 2.0;
      "Se", 2.4;
      "Br", 2.8;
      "Rb", 0.8;
      "Sr", 1.0;
      "Y", 1.2;
      "Zr", 1.4;
      "Nb", 1.6;
      "Mo", 1.8;
      "Tc", 1.9;
      "Ru", 2.2;
      "Rh", 2.2;
      "Pd", 2.2;
      "Ag", 1.9;
      "Cd", 1.7;
      "In", 1.7;
      "Sn", 1.8;
      "Sb", 1.9;
      "Te", 2.1;
      "I", 2.5;
      "Cs", 0.7;
      "Ba", 0.9;
      "La", 1.0;
      (* Skipping the elements between La and Lu *)
      "Hf", 1.3;
      "Ta", 1.5;
      "W", 1.7;
      "Re", 1.9;
      "Os", 2.2;
      "Ir", 2.2;
      "Pt", 2.2;
      "Au", 2.4;
      "Hg", 1.9;
      "Tl", 1.8;
      "Pb", 1.9;
      "Bi", 1.9;
      "Po", 2.0;
      "At", 2.2;
      "Fr", 0.7;
      "Ra", 0.9;
      "Ac", 1.1;
      "Th", 1.3;
      "Pa", 1.4;
      "U", 1.4;
      "Np", 1.3 (* skipping the elements between Np and No *);
    ]
    |> String.Map.of_alist_exn
end

(**
  Note: to calculate the total amount of energy stored in a
  molecule. Sum the bond energies for all of the molecule bounds.
  To calculate the amount of heat released/absored to form a molecule.
  
  1. Calculate the amount of energy in the reagent molecules (the
    energy needed to break these molecules into individual atoms)
  2. Calculate the amount of energy in the output molecule
  
  Subtract 1 from 2.
*)
module Heat_of_formation = struct
  (**
    The energy (measured in kJ mole^-1) of single convalent bonds
    between various atoms and theirselves. For instance H-H.

    See: Pauling, Linus. "General Chemistry." 1970. Appendix VIII.
    Table VIII-1. 
  *)
  let single_bond_energies =
    [
      "H", 436.;
      "Li", 111.;
      "Na", 75.;
      "K", 55.;
      "Rb", 52.;
      "Cs", 45.;
      "B", 225.;
      "C", 344.;
      "Si", 187.;
      "Ge", 157.;
      "Sn", 143.;
      "N", 159.;
      "P", 217.;
      "As", 134.;
      "Sb", 126.;
      "Bi", 105.;
      "O", 143.;
      "S", 266.;
      "Se", 184.;
      "Te", 168.;
      "F", 158.;
      "Cl", 243.;
      "Br", 193.;
      "I", 151.;
    ]
    |> String.Map.of_alist_exn

  (**
    Given a covalent bond of two atoms A and B, this function
    returns the energy (kJ mole^-1) of a bond between A and B.

    Accepts four arguments: [energy0], the energy (Kj mole^-1)
    of a single covalent A to A bond; [energy1], the energy of a
    single covalent B to B bond; [x0], the electronegativity of an
    A molecule; and [x1], the electronegativity of a B molecule.
  
    Note: this function returns an empirical approximation that is
    generally accurate to within 3 kJ mole^-1.

    WARNING: when breaking bonds, generally, the amount of energy
    needed will differ. For example: to break the H-O bond in
    hydrogen per oxide, takes 502 kJ mole^-1. Whereas breaking the
    first H-O bond in water takes 424 kJ mole^-1.
  *)
  let bond_energy ~energy0 ~energy1 ~x0 ~x1 =
    ((energy0 + energy1) / 2.0) + (100.0 * square (x0 - x1)) - (6.5 * int_pow (x0 - x1) 4)

  let%expect_test "bond_energy" =
    bond_energy ~energy0:193.0 ~energy1:143.0 ~x0:2.8 ~x1:3.5 |> printf "%0.4f";
    [%expect {| 215.4394 |}]

  (**
    Accepts two strings [x] and [y] that represent the names of two
    elements, such as "O" and "H", and returns an approximation of
    the amount of energy (measured in kJ mole^-1) stored within a
    single covalent bond between them.

    Note that this approximation is generally accurate to within
    approximately 3%.

    Note also that double and triple bonds are generally stronger
    than would be calculated by doubling the value returned here
    (sometimes significantly less).
  *)
  let bond_energy_simple x y =
    let open Electronegativity in
    bond_energy
      ~energy0:(Core.Map.find_exn single_bond_energies x)
      ~energy1:(Core.Map.find_exn single_bond_energies y)
      ~x0:(Core.Map.find_exn tbl x) ~x1:(Core.Map.find_exn tbl y)

  let%expect_test "bond_energy_simple" =
    bond_energy_simple "O" "H" |> printf "%0.4f";
    [%expect {| 460.5296 |}]

  (**
    Accepts one argument: [bonds], a list of the bonds within a
    molecule; and returns an estimate for the total amount of energy
    (measured in kJ/mole) stored within molecule's bonds.

    Note: this should return a first-order approximation of the
    energy stored within a molecule. Generally, this function will
    be accurate to within 10%. To increase accuracy, you will need
    to account for resonance structures (molecules fluctuate between
    bond structures), orbital spin alignments, etc.

    WARNING: this function should only be used for molecules where
    each pairwise bond is a single covalent bond.

    Example: Br_2 O has the following structure:

      Br - O - Br

    [molecule_enthalpy [("Br", "O"); ("Br", "O")]]
  *)
  let molecule_enthalpy = List.sum (module Float) ~f:(Tuple2.uncurry bond_energy_simple)

  let%expect_test "molecule_enthalpy" =
    molecule_enthalpy [ "Br", "O"; "Br", "O" ] |> printf "%0.4f";
    [%expect {| 430.8787 |}]

  let is_unstable qf = qf < 0.0
end
