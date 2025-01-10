open! Core

(** An intuitive version of the string slice function. *)
let slice s start = function
| 0 -> ""
| len -> String.slice s start len

let get_index ~equal x =
  Fn.compose (Option.map ~f:fst) @@ Array.findi ~f:(fun _ -> equal x)

(**
  Accepts four arguments:
  * [f], a function that accepts two arguments [acc] and [i],
    an integer, and returns either Continue of Stop with a value
  * [finish], a function that accepts a value and returns the
    final value
  * [init], an initial value
  * and [n], an integer
  and folds [f] over the sequence of numbers 0 .. [n] - 1 starting
  with [init] and finalized using [finish].
*)
let loop_until ~f ~finish init n = Sequence.init n ~f:Fn.id |> Sequence.fold_until ~init ~f ~finish

(**
  Accepts three arguments:
  * [f], a function that accepts two arguments [acc] and [i],
    an integer, and returns either Continue or Stop with a value
  * [init], an initial value
  * and [n], an integer
  and performs like a for loop that accumulates [f] over the
  numbers 0 .. [n] - 1 starting with [init] and stopping early iff [f]
  returns None.
*)
let loop = loop_until ~finish:Fn.id

let%expect_test "loop" =
  loop ~f:(fun x y -> Continue (x + y)) 0 5 |> printf !"%d";
  [%expect {| 10 |}]