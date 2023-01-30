# BellBruno

<!--[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://stephans3.github.io/BellBruno.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://stephans3.github.io/BellBruno.jl/dev/)
[![Build Status](https://github.com/stephans3/BellBruno.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/stephans3/BellBruno.jl/actions/workflows/CI.yml?query=branch%3Amain)-->

Compute [Bell polynomials](https://en.wikipedia.org/wiki/Bell_polynomials) for [Faà di Bruno's formula](https://en.wikipedia.org/wiki/Fa%C3%A0_di_Bruno's_formula).

**Basic usage**

```julia
N_der = 10;              # Maximum order of Bell polynomial
bp = bell_poly(N_der);   # Create bell polynomials
bc = bell_coeff(bp);     # Compute bell coefficients
```


**Computing and saving Bell polynomials**

If you want to compute and save the Bell polynomials you may use

```julia
N_der = 20;
path_2_folder = "my_folder_path"
bp = bell_poly( N_der; 
                save_on_disk=true, 
                path_to_folder=path_2_folder, 
                print_iteration=true)
```

or simply
```julia
bp = bell_poly(N_der, save_on_disk = true)
```
to save the files in folder `"bell_results/"`.


**Reading the Bell polynomials**

The Bell polynomials can be read with

```julia
bp_new = read_bell_poly(path_to_folder=path_2_folder)
```

or if the standard folder `"bell_results/"` is used

```julia
bp_new = read_bell_poly()
```


