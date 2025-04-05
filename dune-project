(lang dune 3.16)
(name chemistry)
(license GPL-3.0)
(authors "Larry D. Lee Jr. <llee454@gmail.com>")
(maintainers "Larry D. Lee Jr. <llee454@gmail.com>")
(source (github llee454/chemistry))
(version 0.0.1)
(generate_opam_files)
(package
 (name chemistry)
 (synopsis "A library for performing chemistry calculations")
 (description "A library for performing chemistry calculations based on Linus Pauling's \"General Chemistry\".")
 (depends
  core
  dune
  ocaml-lsp-server
  ocamlformat
  ppx_jane
  ppx_here
  eio_main
  (ocaml_math (= 2.1.0))
  conf-gsl
 )
)
(using coq 0.2)
