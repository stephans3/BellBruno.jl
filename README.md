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
