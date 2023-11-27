# BellBruno

<!--[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://stephans3.github.io/BellBruno.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://stephans3.github.io/BellBruno.jl/dev/)
[![Build Status](https://github.com/stephans3/BellBruno.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/stephans3/BellBruno.jl/actions/workflows/CI.yml?query=branch%3Amain)-->

Compute [Bell polynomials](https://en.wikipedia.org/wiki/Bell_polynomials) for [Faà di Bruno's formula](https://en.wikipedia.org/wiki/Fa%C3%A0_di_Bruno's_formula).

## Basic usage

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

## Quick Tutorial

We are wish to find the derivatives of function composition

$$ h(t) = f(g(t)) = \exp(p~\sin(t))$$

up to the order of $N=10$. We assume $p=0.1$. 

The derivative of the outer function is noted as

$$ \frac{d^n}{dx^n} f(x) = p^{n} \exp(p~x)$$

and the derivative of the inner function is noted as

$$ \frac{d^n}{dt^n} g(t) = \sin(t+n~\pi/2).$$

**Derivative of outer function**

```julia
p = 0.1;
f_der(x,n) = p^n * exp(p*x) 
```

**Inner function and its derivatives**

```julia
g(t) = sin(t)
g_der(t,n) = sin(t+n*π/2)
```

**Bell polynomial and coefficients data**

```julia
using BellBruno
bp = bell_poly(10)
bc = bell_coeff(bp)
```

**Sampling points and Faà di Bruno's formula**
```julia
tgrid = -π:0.01:2π # Sampling points
diff_data = faa_di_bruno(f_der, g, g_der, tgrid,bp, bc)
```

The full listing [tutorial_example.jl can be found here](https://github.com/stephans3/BellBruno.jl/blob/master/examples/src/tutorial_example.jl).

**Derivatives up to order 5**

<img src="https://raw.githubusercontent.com/stephans3/BellBruno.jl/main/examples/results/images/tutorial_1.png" width="450" height="300">

**Derivatives up to order 10**

<img src="https://raw.githubusercontent.com/stephans3/BellBruno.jl/main/examples/results/images/tutorial_2.png" width="450" height="300">