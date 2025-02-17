# MT.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ayushinav.github.io/MT.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ayushinav.github.io/MT.jl/dev/)
[![Build Status](https://travis-ci.com/ayushinav/MT.jl.svg?branch=main)](https://travis-ci.com/ayushinav/MT.jl)
[![Coverage](https://codecov.io/gh/ayushinav/MT.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/ayushinav/MT.jl)
[![Coverage](https://coveralls.io/repos/github/ayushinav/MT.jl/badge.svg?branch=main)](https://coveralls.io/github/ayushinav/MT.jl?branch=main)

`MT.jl` is supposed to be a high performance code for doing forward and inverse modeling in geophysics using julia. We hope to write the code structure such that any other geophysical survey can also be used and we can tend towards a joint forward and inverse modeling library.

## Forward modeling

While forward modeling typically requires solving a PDE obtained using the quasi-static approximation, in 1D, we are fortunate to have the solution for surface impedance in a more analytical form. Currently, this is what is supported.

Supported methods:
* 1D Magnetotellurics (MT)

## Inverse modeling
No surprises here that we are almost always trying to solve for an under-determined system.

Deterministic schemes supported:
* Occam
* Nonlinear schemes using NonlinearSolve.jl 
* Nonlinear schemes using Optimization.jl

Probabilistic schemes supported:'
* MCMC with fixed grids
* MCMC with flexible grids
* RTO-TKO

## Rock physics
* [SEO3](@ref SEO3) : Constable, 2006 : [SEO3: A new model of olivine electrical conductivity](https://doi.org/10.1111/j.1365-246X.2006.03041.x)
* [UHO2014](@ref UHO2014) : Gardés et al., 2014 : [Toward a unified hydrous olivine electrical conductivity law](https://doi.org/10.1002/2014GC005496)
* [Jones2012](@ref Jones2012) : Jones et al., 2012 : [Calibrating laboratory-determined models of electrical conductivity of mantle minerals using geophysical and petrological observations](https://doi.org/10.1029/2012GC004055)
* [Poe2010](@ref Poe2010) : Poe et al., 2010 : [Electrical conductivity anisotropy of dry and hydrous olivine at 8GPa](https://doi.org/10.1016/j.pepi.2010.05.003)
* [Wang2006](@ref Wang2006) : Wang et al., 2006 : [The effect of water on the electrical conductivity of olivine](https://doi.org/10.1038/nature05256)
* [Yoshino2009](@ref Yoshino2009) : Yoshino et al., 2009 :[The effect of water on the electrical conductivity of olivine aggregates and its implications for the electrical structure of the upper mantle](https://doi.org/10.1016/j.epsl.2009.09.032)
* [const_matrix](@ref const_matrix) : The solid matrix is assumed to be of constant bulk conductivity

### Melt
* [Ni2011](@ref Ni2011) : Ni et al., 2011 : [Electrical conductivity of hydrous basaltic melts: implications for partial melting in the upper mantle](https://doi.org/10.1007/s00410-011-0617-4)
* [Sifre2014](@ref Sifre2014) : Sifre et al., 2014 : [Electrical conductivity during incipient melting in the oceanic low-velocity zone](https://doi.org/10.1038/nature13245)
* [Gaillard2008](@ref Gaillard2008) : Gaillard et al., 2008 : [Carbonatite Melts and Electrical Conductivity in the Asthenosphere](https://www.science.org/doi/10.1126/science.1164446)
