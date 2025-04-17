# Rock physics

## Models

We support the following rock physics models to obtain the electrical conductivity of the rock:

### Minerals

  - [SEO3](@ref SEO3) : Constable, 2006 : [SEO3: A new model of olivine electrical conductivity](https://doi.org/10.1111/j.1365-246X.2006.03041.x)
  - [UHO2014](@ref UHO2014) : Gardés et al., 2014 : [Toward a unified hydrous olivine electrical conductivity law](https://doi.org/10.1002/2014GC005496)
  - [Jones2012](@ref Jones2012) : Jones et al., 2012 : [Calibrating laboratory-determined models of electrical conductivity of mantle minerals using geophysical and petrological observations](https://doi.org/10.1029/2012GC004055)
  - [Poe2010](@ref Poe2010) : Poe et al., 2010 : [Electrical conductivity anisotropy of dry and hydrous olivine at 8GPa](https://doi.org/10.1016/j.pepi.2010.05.003)
  - [Wang2006](@ref Wang2006) : Wang et al., 2006 : [The effect of water on the electrical conductivity of olivine](https://doi.org/10.1038/nature05256)
  - [Yoshino2009](@ref Yoshino2009) : Yoshino et al., 2009 :[The effect of water on the electrical conductivity of olivine aggregates and its implications for the electrical structure of the upper mantle](https://doi.org/10.1016/j.epsl.2009.09.032)
  - [const_matrix](@ref const_matrix) : The solid matrix is assumed to be of constant bulk conductivity

### Melt

  - [Ni2011](@ref Ni2011) : Ni et al., 2011 : [Electrical conductivity of hydrous basaltic melts: implications for partial melting in the upper mantle](https://doi.org/10.1007/s00410-011-0617-4)
  - [Sifre2014](@ref Sifre2014) : Sifre et al., 2014 : [Electrical conductivity during incipient melting in the oceanic low-velocity zone](https://doi.org/10.1038/nature13245)
  - [Gaillard2008](@ref Gaillard2008) : Gaillard et al., 2008 : [Carbonatite Melts and Electrical Conductivity in the Asthenosphere](https://www.science.org/doi/10.1126/science.1164446)

## Forward calculations

Forward calculations are fairly easy and involve calling `forward` function on the model. Following is an example on how we can get the conductivity using `Poe2010`:

```julia
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

```@example rp

mix1 = construct_mixing_models([1000 + 273.0, 2e4], [:T, :Ch2o_m],
    [0.1], [SEO3, Ni2011], [HS1962_plus()])

log_cond_mix = forward(mix1, [])
```

**Note**:

!!! note
    

Even when you have a single phase, the use of [`construct_mixing_models`](@ref construct_mixing_models) is recommended using `single_phase` mixing scheme. Make sure to then have `ϕ = 1.`.

```@example rp

mix_single = construct_mixing_models([1000 + 273.0, 2e4], [:T, :Ch2o_ol],
    [1.0], [UHO2014], [single_phase()])

log_cond = forward(mix_single, [])
```

## Stochastic inversion

Now that we know how to compute the rock physics responses, we move towards performing stochastic inversion for rock physics parameters given some rock conductivity and the corresponding uncertainty.

```julia
using Distributions
using StatsPlots
using Turing
# define a rock physics distribution to sample parameters from
m_dist = RockphyModelDistribution(
    Product([
        Uniform(1173.0, 1573.0),
        Uniform(1e2, 3e4)
    ]),
    [:T, :Ch2o_m],
    Product([Uniform(0.0, 0.4)]),
    [SEO3, Ni2011],
    [HS1962_plus()]
)
```

[`RockphyModelDistribution`](@ref RockphyModelDistribution) allows us to define the a priori distribution for rock physics parameters. If the user wants to keep a certain parameter constant, they should define the range of *a priori* distribution as very small, eg., if, in the above, we want to keep water content in melt `CH2o_m` constant, we define it's *a apriori* as a narrow range as

```julia
m_dist = RockphyModelDistribution(
    Product([
        Uniform(1173.0, 1573.0),
        Uniform(2e4 - 1, 2e4 + 1) # Keeping Ch2o_m almost fix around 2e4
    ]),
    [:T, :Ch2o_m],
    Product([Uniform(0.0, 0.4)]),
    [SEO3, Ni2011],
    [HS1962_plus()]
)
```

```@example rp
r_dist = RockphyResponseDistribution((x, y) -> MultivariateNormal(x, y))

m_cache = mcmc_cache(
    m_dist,
    r_dist,
    50_000,
    MH()
)

# rock physics response
rp_obs = RockphyCond(log_cond_mix.σ)
err_rp = RockphyCond(log_cond_mix.σ .* 0.01)

samples = stochastic_inverse(rp_obs, err_rp, [], m_cache) #, trans_utils=(m=log_tf,))

plot(samples)
```
