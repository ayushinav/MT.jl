"""
    two_phase_modelDistributionType(m1, m2, mix)

Rock physics model distribution type to combine two phases.

## Arguments

  - `m1` : model distribution type corresponding to phase 1
  - `m2` : model distribution type corresponding to phase 2
  - `mix` : mixing type, available options are `HS_1962_plus()`, `HS1962_minus`, `MAL(m)`

## Usage

```julia
two_phase_modelType(SEO3Distribution, Ni2011Distribution, HS1962_plus())
```
"""
mutable struct two_phase_modelDistributionType{T1, T2, M}
    m1::Type{T1}
    m2::Type{T2}
    mix::M
end

"""
    two_phase_modelDistribution(ϕ, m1, m2, mix)

Rock physics model distribution to combine two phases, usually constructed through `two_phase_modelDistributionType`[@ref]

## Arguments

  - `ϕ` : distribution (or value) of vol. fraction of the **second** phase
  - `m1` : model distribution corresponding to phase 1
  - `m2` : model distribution corresponding to phase 2
  - `mix` : mixing type, available options are `HS_1962_plus()`, `HS1962_minus`, `MAL(m)`

## Usage

```julia
m = two_phase_modelType(SEO3Distribution, Ni2011Distribution, HS1962_plus())
ps_nt_dist = (; T=product_distribution(Uniform(1200.0f0, 1400.0f0)),
    Ch2o_m=MvNormal([100.0f0], diagm([20.0f0])), ϕ=[0.1f0])
model = m(ps_nt_dist)

resp = forward(model)
```
"""
mutable struct two_phase_modelDistribution{V, T1, T2, M} <: AbstractRockphyModelDistribution
    ϕ::V
    m1::T1
    m2::T2
    mix::M
end

two_phase_modelDistributionType(m1) = m1
two_phase_modelDistributionType(m1, m::phase_mixing) = m1

function (model::two_phase_modelDistributionType)(ps::NamedTuple)
    mix = model.mix
    ϕ = ps.ϕ

    v1 = from_nt(model.m1, ps)

    v2 = from_nt(model.m2, ps)

    return two_phase_modelDistribution(ϕ, v1, v2, mix)
end

function from_nt(m::Type{T}, nt::NamedTuple) where {T <: two_phase_modelDistributionType}

    # @show m.types[2].parameters[1]
    # @show m.types[3] #[3].parameters[1]
    ϕ = nt.ϕ
    m1 = m.types[1].parameters[1]
    m2 = m.types[2].parameters[1]
    mix = m.types[3] #.parameters[1]

    model1 = MT.from_nt(m1, nt)
    model2 = MT.from_nt(m2, nt)

    return two_phase_modelDistribution(ϕ, model1, model2, mix())
end