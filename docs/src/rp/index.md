# Rock physics

## Models

We support the following rock physics models to obtain the electrical conductivity of the rock:

### Minerals
[SEO3: A new model of olivine electrical conductivity](https://doi.org/10.1111/j.1365-246X.2006.03041.x)

  - [UHO2014](@ref UHO2014) : [Toward a unified hydrous olivine electrical conductivity law](https://doi.org/10.1002/2014GC005496)
  - [Jones2012](@ref Jones2012) : [Calibrating laboratory-determined models of electrical conductivity of mantle minerals using geophysical and petrological observations](https://doi.org/10.1029/2012GC004055)
  - [Poe2010](@ref Poe2010) : 
  - [Wang2006](@ref Wang2006) : [The effect of water on the electrical conductivity of olivine](https://doi.org/10.1038/nature05256)
  - [Yoshino2009](@ref Yoshino2009) :[The effect of water on the electrical conductivity of olivine aggregates and its implications for the electrical structure of the upper mantle](https://doi.org/10.1016/j.epsl.2009.09.032)
  - [const_matrix](@ref const_matrix) : The solid matrix is assumed to be of constant bulk conductivity

### Melt

  - [Ni2011](@ref Ni2011) : [Electrical conductivity of hydrous basaltic melts: implications for partial melting in the upper mantle](https://doi.org/10.1007/s00410-011-0617-4)
  - [Sifre2014](@ref Sifre2014) : [Electrical conductivity during incipient melting in the oceanic low-velocity zone](https://doi.org/10.1038/nature13245)
  - [Gaillard2008](@ref Gaillard2008) : [Carbonatite Melts and Electrical Conductivity in the Asthenosphere](https://www.science.org/doi/10.1126/science.1164446)

## Forward calculations

Forward calculations are fairly easy and involve calling `forward` function on the model. Following is an example on how we can get the conductivity using `Poe2010`:

```@example rp
using MT
model_poe = Poe2010(1000 + 273.0, 2e4)

log_cond = forward(model_poe, [])
```

### Mixing phases

We can also mix in melt and use the following mixing schemes to get the bulk conductivity:

  - [HS1962_plus](@ref HS1962_plus) : Hashim-Strikman upper bound for 2 phases
  - [HS1962_minus](@ref HS1962_minus) : Hashim-Strikman lower bound for 2 phases
  - [MAL](@ref MAL) : Modified Archie's law
  - [single_phase](@ref single_phase) : Bulk rock is composed of one material

Below we estimate the bulk conductivity of rock with solid phase conductivity governed by `SEO3` and melt by `Ni2011`. We'll assume a porosity of 0.1 and calculate the upper HS bounds for the matrix.

```julia

mix1 = construct_mixing_models([1000 + 273.0, 2e4], [:T, :Ch2o_m],
    [0.1], [SEO3, Ni2011], [HS1962_plus()])

log_cond_mix = forward(mix1, [])
```

**Note**:

!!! note
    

Even when you have a single phase, the use of [`construct_mixing_models`] **<===  NEED REF HERE** is recommended using `single_phase` mixing scheme. Make sure to then have `ϕ = 1.`.

```julia

mix_single = construct_mixing_models([1000 + 273.0, 2e4], [:T, :Ch2o_ol],
    [1.0], [UHO2014], [single_phase()])

log_cond = forward(mix_single, [])
```
