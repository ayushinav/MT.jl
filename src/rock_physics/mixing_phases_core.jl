abstract type phase_mixing end

"""
    HS1962_plus

Hashim-Strikman upper bound for mixing 2 phases

## Usage

```julia
model_mix = construct_mixing_models([1000.0 + 273.0, 2e4], [:T, :Ch2o_m],
    Product([Uniform(0.0, 1.0)]), [SEO3, Ni2011], [HS1962_plus()])
log_cond = forward(model_mix, [])
```

## References

  - Paul W. J. Glover (2010), "A generalized Archie's law for n phases",
    Geophysics 2010; 75 (6): E247–E265. doi: https://doi.org/10.1190/1.3509781
"""
mutable struct HS1962_plus <: phase_mixing end

"""
    HS1962_minus

Hashim-Strikman lower bound for mixing 2 phases

## Usage

```julia
model_mix = construct_mixing_models(
    [1000.0 + 273.0, 2e4], [:T, :Ch2o_m], [0.1], [SEO3, Ni2011], [HS1962_minus()])
log_cond = forward(model_mix, [])
```

## References

  - Paul W. J. Glover (2010), "A generalized Archie's law for n phases",
    Geophysics 2010; 75 (6): E247–E265, doi: https://doi.org/10.1190/1.3509781
"""
mutable struct HS1962_minus <: phase_mixing end

"""
    MAL

Modified Archie's law for mixing 2 phases

## Usage

```julia
model_mix = construct_mixing_models([1000. + 273., 2e4]
    [:T, :Ch2o_m]],
    Product([Uniform(0., 1.)]),
    [SEO3, Ni2011],
    [MAL(1.2)]
)
log_cond = forward(model_mix, [])
```

## References

- Glover, P. W. J., Hole, M. J., & Pous, J. (2000),
"A Modified Archie’s Law for two conducting phases",
Earth and Planetary Science Letters, 180(3–4), 369–383, doi: https://doi.org/10.1016/S0012-821X(00)00168-0
"""
mutable struct MAL <: phase_mixing
    m
end

mutable struct HSn_plus <: phase_mixing end

mutable struct HSn_minus <: phase_mixing end

"""
    single_phase

Single phase only conductivity. Assumes the rock matrix is composed of a single phase only.
"""
mutable struct single_phase <: phase_mixing end

# mixing functions

function mix_models(σs, ϕ, ::HS1962_plus)
    σ_max = maximum(σs)
    σ_min = minimum(σs)
    phi = first(ϕ)

    num = 3 * (1 - phi) * (σ_max - σ_min) # numerator
    den = 3 * σ_max - phi * (σ_max - σ_min) # denominator
    esig = σ_max * (1 - (num / den))

    return esig
end

function mix_models(σs, ϕ, ::HS1962_minus)
    σ_max = maximum(σs)
    σ_min = minimum(σs)
    phi = first(ϕ)

    num = 3 * (phi) * (σ_max - σ_min) # numerator
    den = 3 * σ_min + (1 - phi) * (σ_max - σ_min) # denominator
    esig = σ_min * (1 + (num / den))

    return esig
end

function mix_models(σs, ϕ, mal::MAL)
    σ_fluid = (σs[2])
    σ_matrix = (σs[1])

    phi = first(ϕ)
    sig = σ_fluid

    if phi < 1
        p = log10(1 - phi^mal.m) * inv(log10(1 - phi))
        sig = σ_fluid * phi^mal.m + σ_matrix * (1 - phi)^p
    end

    return sig
end

function mix_models(σs, ϕ, ::single_phase)
    @assert length(σs) == 1
    return first(σs)
end
