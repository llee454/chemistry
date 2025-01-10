Chemistry Readme
=========================

This library defines a collection of functions for performing chemistry calculations. The primary reference for this library is Linus Pauling's classic text book "General Chemistry." 

The OCaml and Coq code use Dune and OPAM. You can configure and compile these files using the following comamnds:

```bash
opam switch create . 5.1.0+options --no-install
eval $(opam env)
opam update
opam upgrade
dune build chemistry.opam # to generate OPAM package file
opam install --deps-only . -y
dune build
dune runtest
dune exec src/main.exe
```