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

```julia
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
