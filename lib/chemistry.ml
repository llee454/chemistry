(**
  The equations presented here are from:

  1. Pauling, Linus. 1988. "General Chemistry." Dover.
*)
open! Core
open Float
open! Ocaml_math
 
include Coq_chemistry
 
(** The number of meters spanned by an angstrom *)
let angstrom = 1e-10

let angstroms_to_meters = ( * ) angstrom

let meters_to_angstroms x = x / angstrom

(** Planck's constant measured in Joule second *)
let plancks_constant = 6.6252e-34

let planck_bar = plancks_constant / (2. * pi)

(** The speed of light measured in meters per second. *)
let light_velocity = 299_792_458.

(** The avogrado's number of atoms - number of atoms in a mole. *)
let mole = 6.0229e23

(** The mass of a proton in kilograms. *)
let proton_mass = 1.67239e-27

(** The mass of a neutron in kilograms. *)
let neutron_mass = 1.67470e-27

(** The mass of a dalton in kilograms. *)
let dalton_to_kg = 1.66033e-27

let celcius_to_kelvin x = x + 273.15

let kelvin_to_celcius x = x - 273.15

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
let electron_mass_kg = 9.1083e-31

(** The charge of an electron in coulombs. *)
let electron_charge_c = 1.60206e-19

let electron_charge_s = electron_charge_c * columnb_to_stoney

(** The number of joules represented by one electron volt *)
let electron_volt = 1.60206e-19

let electron_volts_to_joules = ( * ) electron_volt

let joules_to_electron_volts x = x / electron_volt

(**
  Note that j L^-1 = kN m^-2
*)
let atmospheres_to_joule_liters atm = 101.325 * atm

(**
  Accepts the [mass], measured in Kg, and [energy], measured in j, of a
  moving particle; and returns its velocity in m/s.
*)
let get_velocity ~mass ~energy =
  sqrt (2.0 * energy/mass)

(**
  Accepts the [mass], measured in Kg, and [velocity], measured in m/s,
  of a particle and returns the energy of that particle measured in joules.
*)
let get_energy ~mass ~velocity =
  mass * (square velocity)/2.0

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

    WARNING: in [1] Pauling gives a value of 2.02 Å for the unit test
    on page 143. This function may be incorrect, but I cannot find any
    discrepancy between it and the equation printed in the book. I assume
    that the text itself contains a misprint.
  *)
  let closest_radius ?electron_mass ~z ~n l =
    (square n)/z * (bohr_radius_constant ?electron_mass ()) * (1.5 - ((square l + l)/(2. * square n)))

  let%expect_test "closest_radius" =
    (closest_radius ~electron_mass:(reduced_electron_mass (19.0 * proton_mass)) ~z:1. ~n:3. 2.) / angstrom |> printf !"%0.4f Å\n";
    [%expect {| 5.5565 Å |}]

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

  module Ionic = struct

    (**
      This array represents the empirically observed correlation between the
      electronegativity difference between two atoms linked by a covalent bond
      and the degree to which that bond is ionic.

      See [1] Table 6.5 page 184.
    *)
    let empirical_degree_ionic =
      [|
        [| 0.2;  0.01 |];
        [| 0.4;  0.04 |];
        [| 0.6;  0.09 |];
        [| 0.8;  0.15 |];
        [| 1.0;  0.22 |];
        [| 1.2;  0.3  |];
        [| 1.4;  0.39 |];
        [| 1.6;  0.47 |];
        [| 1.8;  0.55 |];
        [| 2.0;  0.63 |];
        [| 2.2;  0.7  |];
        [| 2.4;  0.76 |];
        [| 2.6;  0.82 |];
        [| 2.8;  0.86 |];
        [| 3.0;  0.89 |];
        [| 3.2;  0.92 |];
      |]


    let get_coeffs =
      let n = Array.length empirical_degree_ionic
      and f Nonlinear_fit.{ks; x} =
        match ks with
        | [| k0; k1 |] ->
          k0 * cdf_gaussian_p ~x:(log (k1 * x)) ~std:1.0
        | _ -> failwiths ~here:[%here] "Error: an internal error occured." ks [%sexp_of: float array]
      in
      let xs = Array.init n ~f:(fun i -> empirical_degree_ionic.(i).(0))
      and ys = Array.init n ~f:(fun i -> empirical_degree_ionic.(i).(1))
      in
      Nonlinear_fit.f ~f ~ks_init:[|1.0; 1.0|] ~xs ~ys 

    let%expect_test "get_coeffs" =
      printf !"%{sexp: float array}" get_coeffs;
      [%expect {| (1.8615176504138462 0.32202559402545522) |}]

    (**
      Accepts one argument: [electronegativity_diff] which represents the
      absolute value of the difference in eletronegativity between two atoms
      linked by a covalent bond; and returns an approximation of the degree
      to which that bond is ionic - i.e. the degree to which the shared
      electrons disproportionately orbit the more electronegative of the
      two atoms.contents

      Note: this function is a regression fit for the empirical table presented
      in [1] Table 6.5 page 184. It is accurate to within 3.4% for all values
      less than and equal to 3.2 and the error increases as the difference
      in electronegativity increases.
    *)
    let get_approx_degree_ionic electronegativity_diff =
      if [%equal: float] electronegativity_diff 0.0
      then 0.0
      else
        let k0 = 1.8615176504138462
        and k1 = 0.32202559402545522
        in
        k0 * cdf_gaussian_p ~x:(log (k1 * electronegativity_diff)) ~std:1.0

    let%expect_test "get_approx_degree_ionic" =
      Array.map empirical_degree_ionic ~f:(function
        | [|electronegativity_diff; degree_ionic|] ->
          abs (degree_ionic - get_approx_degree_ionic electronegativity_diff)
        | xs -> failwiths ~here:[%here] "Error: an internal error occured." xs [%sexp_of: float array])
      |> printf !"%{sexp: float array}";
      [%expect {|
        (0.0043259164570501416 0.0023772825511848006 0.0032495216360882706
         0.012895802386018951 0.019355873781161714 0.018044520067233949
         0.0061804344766568375 0.0021304889854154019 0.0050208292778404218
         0.015742256620153605 0.020225121515539013 0.018489778073564778
         0.020452977300009612 0.0059707800194609417 0.015134013368263921
         0.03305434060395529)
        |}]
  end
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

(**
  This module defines functions that use the electroneutrality principle to
  devise chemical formula for molecules.

  Given a set of n atoms that comprise a molecule [a_0, a_1, ..., a_{n-1}],
  we can:
  1. assign unshared electrons to individual atoms
  2. place zero or more covalent bonds (shared pairs of electrons) between two atoms.

  Every atom is described by two properties:
  1. the number of valence electrons (the number of electrons in its outermost/"valence" shell)
  2. its electronegativity

  We use a symmetric matrix to count the number of shared electrons between each pair of atoms.

  [
    [y_{0,0}, y_{0,1}, ..., y_{0,n-1}],
    [y_{1,0}, y_{1,1}, ..., y_{1,n-1}],
    ...
    [y_{n-1,0}, y_{n-1,1}, ..., y_{n-1,n-1}],
  ]

  The matrix is symmetric so the y_{i,j} = y_{j,i} for all i and j. Also, we assume that all diagonal elements y_{i,i} represent the number of unshared electrons that we assign to the ith atom.

  Because the matrix is symmetric, we only consider the top right triangular entries. These can be listed in a single vector as:

  [y_{0,0}, y_{0,1}, ..., y_{0,n-1}, y_{1,1}, ... y{1,n-1}, ..., y{n-1,n-1}]

  We call this vector the "electron configuration vector".

  Given an "electron configuration vector," we can compute the charge on each atom. For each atom a_i, we calculate m_i, the number of unshared electrons and covalent bonds involving a_i and then compute:

    sum (ionic (a_i.electronegativity - a_j.electronegativity), j, n-1) - (m - a_i.num_valence_electrons)

  where ionic (x, y) returns the degree to which an electron will favor an electron with electronegativity x over y. 

  We use a genetic algorithm to find an electron configuration that minimizes the charge across atoms.
*)
module Electroneutrality = struct
  module Atom = struct
    type t = {
      num_valence_electrons: int;
      electronegativity: float; 
    } [@@deriving sexp]
  end

  module Electron_configuration_vector = struct
    (**
      Accepts an electron configuration vector and returns the number of atoms represented.

      n = x(x + 1)/2
      0 = x^2 + x - 2n
    *)
    let get_num_atoms x =
      let n = Array.length x in
      Float.iround_nearest_exn @@ (sqrt (8.0*float n + 1.0) - 1.0)/2.0

    (**
      Accepts three arguments: [n], the number of atoms in an electron
      configuration matrix; [i], a row index; [j], a column index; and returns
      the index of the corresponding element within the associated electron
      configuration vector.
    *)
    let get_index n i j =
      let open Int in
      let k = min i j
      and l = max i j
      in
      k*n - (k*(k + 1))/2 + l

    (*
      [| 0.0; 1.0; 2.0;
              3.0; 4.0;
                   5.0 |]
    *)
    let%expect_test "get" =
      [
        get_index 3 0 0;
        get_index 3 1 1;
        get_index 3 2 2;
        get_index 3 1 2;
        get_index 3 2 0;
        get_index 3 1 0;
        get_index 3 1 2;
      ]
      |> printf !"%{sexp: int list}";
      [%expect {| (0 3 5 4 2 1 4) |}]

    (**
      Accepts an electron configuration vector and returns the corresponding
      electron configuration matrix.
    *)
    let get x i j = x.(get_index (get_num_atoms x) i j)

    (**
      Accepts two arguments: [atoms], a list of atom descriptions; and [xs],
      an electron configuration vector; and returns the square of the charge
      of the atoms in the described molecular configuration.
    *)    
    let get_charge (atoms : Atom.t array) (config : int array) =
      Array.foldi atoms ~init:0.0 ~f:(fun i total_square_charge ai ->
        let (ai_charge, ai_num_electrons) = Array.foldi atoms ~init:(0.0, get config i i)
          ~f:(fun j ((ai_charge, ai_num_electrons) as acc) aj ->
            if [%equal: int] i j
            then acc
            else
              let num_shared_electrons = get config i j
              and partial_charge = Electronegativity.Ionic.get_approx_degree_ionic @@ abs @@ ai.electronegativity - aj.electronegativity
              and charge_sign = if ai.electronegativity < aj.electronegativity
                then 1.0
                else -1.0
              in
              (* printf "charge between %d and %d: %f (%f for %d electrons)\n" i j (charge_sign * partial_charge * float num_shared_electrons) partial_charge num_shared_electrons; *)
              (ai_charge + charge_sign * partial_charge * float num_shared_electrons,
               Int.(ai_num_electrons + num_shared_electrons))
        )
        in
        (* printf "charge on %d: %f\n" i @@ ai_charge + float Int.(ai.num_valence_electrons - ai_num_electrons); *)
        total_square_charge + square (ai_charge + float Int.(ai.num_valence_electrons - ai_num_electrons))
      )
    
    let%expect_test "get_charge" =
      let atoms = [|
        (* H *) Atom.{ num_valence_electrons = 1; electronegativity = 2.1 };
        (* C *) Atom.{ num_valence_electrons = 4; electronegativity = 2.5 };
        (* N *) Atom.{ num_valence_electrons = 5; electronegativity = 3.0 }
      |]
      (* H - N --- C: *)
      and config0 = [|
        0; 0; 1;
           2; 3;
              0
      |]
      (* H - C --- N: *)
      and config1 = [|
        0; 1; 0;
           0; 3;
              2
      |]
      in
      [|
        get_charge atoms config0;
        get_charge atoms config1;
      |] |> printf !"%{sexp: float array}";
      [%expect {| (1.0695658090848179 0.060284980182740575) |}]
  end
end


module Gases = struct
  (** The Boltzmann constant "k" measured in j deg^-1. *)
  let boltzmann_constant = 1.3805E-23

  (** Usually denotes "R" measured in j deg^-1 mole^-1 *)
  let gas_constant_joules = boltzmann_constant * mole

  let%expect_test "gas_constant_joules" =
    printf "%f" gas_constant_joules;
    [%expect {| 8.314613 |}]

  (** Usually denoted "R" measured in l atm deg^-1 mole^-1 *)
  let gas_constant_atm_liters = gas_constant_joules / atmospheres_to_joule_liters 1.0

  let%expect_test "gas_constant_atm_liters" =
    printf "%0.04f atm l deg^-1 mole^-1" gas_constant_atm_liters;
    [%expect {| 0.0821 atm l deg^-1 mole^-1|}]

  (** The temperature of the standard atmosphere measured in degrees Kelvin *)
  let standard_temperature = 298.15

  (**
    Accepts two arguments: [num_moles] the number of moles of gas molecules;
    and [temperature], measured in degrees Kelvin; and returns the energy
    measured in joules that an ideal gas pushes outwards.
  *)
  let get_pressure_energy ~num_moles ~temperature =
    num_moles * gas_constant_joules * temperature

  (**
    Accepts three argumnets: [num_moles] the number of moles of the gas
    molecules; [temperature] measured in degrees Kelvin; and [pressure],
    measured in atmospheres; and returns the volume that an ideal gas would
    occupy when subjected to the given pressure.
  *)
  let get_volume ~num_moles ~temperature ~pressure =
    get_pressure_energy ~num_moles ~temperature / pressure

  (**
    Accepts three arguments: [temperature] measured in degrees Kelvin;
    [pressure] measured in atmospheres; [mass] measured in grams per mole
    of the gas molecules; and returns the density of an ideal gas subjected
    to these conditions measured in grams per liter.
  *)
  let get_density ~temperature ~pressure ~mass =
    (pressure * mass)/(gas_constant_atm_liters * temperature)

  (**
    The heat capacity (j deg^-1 mole^-1) of an ideal gas when volume is
    held constant.

    Note: this is typically expressed as 3/2 R.
  *)
  let heat_capacity_volume = 1.5 * gas_constant_joules

  (**
    The heat capacity (j deg^-1 mole^-1) of an ideal gas when the pressure
    is held constant.

    Note: this is typically expressed as 5/2 R.
  *)
  let heat_capacity_pressure = 2.5 * gas_constant_joules

  (**
    Accepts one argument [temperature] (K) and returns the most probable
    kinetic energy (j) of a molecule within an ideal gas having the given
    temperature.

    Note: it is a common mistake to believe that this is kT. Many people,
    including my past self, mistakenly try to calculate this by taking the
    most probable velocity and calculating the corresponding energy, but this
    approach is mistaken without accounting for the rescaling between v and E.
  *)
  let get_most_probable_energy temperature = boltzmann_constant * temperature / 2.0

  let%expect_test "get_most_probable_energy" =
    get_most_probable_energy 400.0 |> joules_to_electron_volts |> printf "%0.4f eV";
    [%expect {| 0.0172 eV |}]

  (**
    Accepts the [mass], measured in Kg, and [temperature], measured in K,
    of a population of molecules, and returns the most probable velocity of
    molecules within the gas, measured in m/s.
  *)
  let get_most_probable_velocity ~mass ~temperature =
    sqrt (2.0 * boltzmann_constant * temperature/mass)

  let%expect_test "get_most_probable_velocity" =
    let helium_molecule_mass_kg = 2.0 * 4.0026/(1_000.0 * mole)
    and temperature = 100.0 in
    get_most_probable_velocity ~mass:helium_molecule_mass_kg ~temperature |> printf "%0.4f";
    [%expect {| 455.7744 |}]

  (**
    Accepts one argument: [temperature] measured in degrees Kelvin; and
    returns the mean kinetic energy `(1/2)mv^2` of the molecules within an
    ideal gas at the given temperature in joules.

    Note: according to [1] page 323, "it has been found that the average
    kinetic energy per molecule, 1/2 mv^2 is the same for all gases at the
    same temperature."
  *)
  let get_mean_energy temperature =
    3.0 * (boltzmann_constant * temperature)/2.0

  let%expect_test "get_mean_energy" =
    let hydrogen_molecule_mass_kg = 2.0 * 1.00797 / (1_000.0 * mole)
    and oxygen_molecule_mass_kg = 2.0 * 15.9994 / (1_000.0 * mole) in
    [
      (hydrogen_molecule_mass_kg, 820.0);
      (oxygen_molecule_mass_kg, 0.0)
    ]
    |> List.map ~f:(fun (mass, temperature) ->
      get_velocity ~mass ~energy:(celcius_to_kelvin temperature |> get_mean_energy)
    )
    |> printf !"%{sexp: float list}";
    [%expect {| (3677.7545172748069 461.44018789203943) |}]

  (**
    Accepts two arguments: [temperature] measured in Kelvins; and [energy]
    measured in electron volts; and returns the probability density that
    a randomly selected molecule from an ideal gas will have the given
    energy.
  *)
  let get_energy_probability ~temperature ~energy =
    let k = 13.805E-24 in
    (*
      boltzmann's constant expressed as electron volts per degree.

      Note: this function works with extremely small numbers which led to
      numerical instability. To improve numerical stability, we switch from
      joules to electron volts.
    *)
    let j = joules_to_electron_volts k in
    let jt = j*temperature in
    (2.0/jt) * sqrt (energy/(pi*jt)) * exp (-energy/jt)

  let%expect_test "get_energy_probability check mean energy" =
    (Integrate.qag () ~f:(fun energy ->
        energy * get_energy_probability ~temperature:standard_temperature ~energy
      ) ~lower:0.0 ~upper:1.00).out
    |> printf "mean energy reference %f eV calculated %f" (get_mean_energy standard_temperature |> joules_to_electron_volts);
    [%expect {| mean energy reference 0.038538 eV calculated 0.038538 |}]

  let%expect_test "get_energy_probability check most probable energy" =
    let temperature = 100.0 in
    let module SA = Simulated_annealing (struct
      type t = float

      let copy x = x

      let energy x = 1.0/get_energy_probability ~temperature ~energy:x

      let step x dist = x +. Random.float (2.0 *. dist) -. dist

      let dist x y = Float.abs (x -. y)

      let print = None
    end)
    in
    let most_probable_energy = SA.(f ~num_iters:1_000 ~step_size:1E-4 (create_state 0.1))
    and most_probable_energy_ref = get_most_probable_energy temperature |> joules_to_electron_volts in
    printf "ref: %f calc %f err %f" most_probable_energy_ref most_probable_energy (abs (most_probable_energy_ref - most_probable_energy));
    [%expect {| ref: 0.004309 calc 0.004309 err 0.000000 |}]

    (**
      Accepts five arguments:
      * [a] the van der Waals' a constant for the gas
      * [b] the van der Waals' b constant for the gas
      * [n] the number of moles
      * [volume] the number of liters
      * [temperature] the temperature (K)
      and uses the van der Waals equation of state to calculate the pressure
      of the gas.
    *)
    let van_der_waals_pressure ~a ~b ~n ~volume ~temperature =
      let r = gas_constant_atm_liters
      in
      (n*r*temperature*(square volume) - a*(square n)*volume + a*b*(pow_int n 3))/
      (pow_int volume 3 - b*n*square volume)

    let%expect_test "van_der_waals_pressure" =
      van_der_waals_pressure ~a:1.39 ~b:0.0391 ~n:1.3 ~volume:2.3 ~temperature:295.13
      |> printf "%f";
      [%expect {| 13.553739 |}]

    (**
      The van_der_waals_volume is cubic in volume.contents

      pv^3 - (nrt+npb)v^2 + n^2av - n^3ab = 0
    *)
    let van_der_waals_volume ~a ~b ~n ~pressure:(p:float) ~temperature:(t:float) =
      let r = gas_constant_atm_liters in
      let result = Polynomial.Cubic.solve
        (- (n*r*t + n*p*b)/p)
        ((square n)*a/p)
        (- (pow_int n 3)*a*b/p)
      in if not @@ [%equal: int] result.n 1
      then failwiths ~here:[%here] "Error: tried to calculate the volume of a gas using the van der Waals equation of state and the equation did not have a unique real solution." (a, b, n, p, t) [%sexp_of: float * float * float * float * float]
      else result.x0

    let%expect_test "van_der_waals_volume" =
      van_der_waals_volume ~a:1.39 ~b:0.0391 ~n:1.3 ~pressure:13.553739 ~temperature:295.13
      |> printf "%f";
      [%expect {| 2.300000 |}]

    let van_der_waals_temperature ~a ~b ~n ~volume ~pressure =
      (pressure + (square n)*a/(square volume))*(volume - n*b)/
      (n*gas_constant_atm_liters)

    let%expect_test "van_der_waals_temperature" =
      van_der_waals_temperature ~a:1.39 ~b:0.0391 ~n:1.3 ~volume:2.3 ~pressure:13.553739
      |> printf "%f";
      [%expect {| 295.130004 |}]
end

module Schrodinger = struct
  (**
    This function evaluates the one-dimensional time-independent shrodinger
    equation (TISE). It accepts five arguments:

    * [mass] (kg) - the mass of the particle 
    * [potential] - a function that accepts a position [x] and returns the
      strength of the field (N) at [x].
    * [energy] (j) - the kinetic energy of the particle
    * [x] (m) - the location of the particle along the x-axis
    * [psi] - the wave function, which accepts a position [x] and returns
      a complex number whose square magnitude represents the probability
      (density) of the particle being at point [x].

    Note if [psi] is a solution of the equation, this function will return 0.

    Note: in the general TISE, [psi] will return a complex number. However,
    this function requires [psi] to return only real numbers.

    Note: the TISE is introduced on page 322, equation 9-14, in [1]
  *)
  let time_independent ~mass ~potential ~energy ~x ~psi =
    -(square plancks_constant)/(8.0*mass*(square pi)) * Deriv.nth ~f:psi ~x ~h:angstrom 2 + (potential x - energy)*(psi x)

  let%expect_test "time_independent box" =
    let n = 10 (* the principal quantum number *)
    and x = 1.0 * angstrom (* the position at which we will evaluate our wave function *)
    and width = 2.0 * angstrom
    and mass = electron_mass_kg
    in
    (* a well potential 2 angstroms wide *)
    let potential x =
      if x <= 0.0 || width <= x
      then 1e6
      else 0.0
    and get_psi n x = sqrt (2.0/width)*sin (pi*(float n)*x/width)
    (* uncomment this to verify that this test will fail when given an invalid wave function. *)
    (* and get_psi _n x = pdf_normal ~mean:(1.0*angstrom) ~std:(0.25*angstrom) x *)
    and get_wavelength n = 2.0*width/(float n)
    and get_energy ~mass ~lambda =
      square (plancks_constant)/(2.0*mass*(square lambda))
    in
    let psi = get_psi n
    and energy = get_energy ~mass ~lambda:(get_wavelength n)
    in
    printf "energy: %f eV, soln: %f"
      (joules_to_electron_volts energy)
      (1e10 * time_independent ~mass ~potential ~energy ~x ~psi);
    [%expect {| energy: 940.008768 eV, soln: -0.000000 |}]

  let%expect_test "time_independent spring" =
    let n = 5 (* the principal quantum number *)
    and x = 10.0 * angstrom (* the position at which we will evaluate our wave function *)
    and k = 1e9 * electron_volt (* the spring constant *)
    and mass = electron_mass_kg
    in
    let w = sqrt (k / mass) (* the angular frequency of the oscillation *)
    and hermite n z =
      match n with
      | 0 -> 1.0
      | 1 -> 2.0*z
      | 2 -> 4.0*(square z) - 2.0
      | 3 -> 8.0*(pow_int z 3) - 12.0*z
      | 4 -> 16.0*(pow_int z 4) - 48.0*(square z) + 12.0
      | 5 -> 32.0*(pow_int z 5) - 160.0*(pow_int z 3) + 120.0*z
      | 6 -> 64.0*(pow_int z 6) - 480.0*(pow_int z 4) + 720.0*(square z) - 120.0
      | _ -> failwiths ~here:[%here] "Error: an error occured while trying to evaluate a hermite polynomial. The requested degree is not unsupported." n [%sexp_of: int]
    in
    (* uncomment this to verify that this test will fail when given an invalid wave function. *)
    (* let get_psi _n x = pdf_normal ~mean:(1.0*angstrom) ~std:(10.0*angstrom) x *)
    (* the solution for the quantum oscillator *)
    let get_psi n x =
      (1.0/sqrt ((pow_int 2.0 n)*(fact n)))
      *(expt (mass*w/(pi*planck_bar)) (1//4))
      *(exp (-mass*w*(square x)/(2.0*planck_bar)))
      *(hermite n (x*sqrt (mass*w/planck_bar)))
    and get_energy n = planck_bar*w*(float n + 0.5)
    in
    let energy = get_energy n
    and psi = get_psi n
    and potential x = -k*x
    in
    printf "energy: %f eV, soln: %f"
      (joules_to_electron_volts energy)
      (1e10 * time_independent ~mass ~potential ~energy ~x ~psi);
    [%expect {| energy: 0.000048 eV, soln: -0.000000 |}]
end