#! format: off

"""
    two_phase_modelType(m1, m2, mix)

Rock physics model to combine two phases.

## Arguments
 - `m1` : model type corresponding to phase 1
 - `m2` : model type corresponding to phase 2
 - `mix` : mixing type, available options are `HS_1962_plus()`, `HS1962_minus`, `MAL(m)`

## Usage

```julia
two_phase_modelType(SEO3, Ni2011, HS1962_plus())
```
"""
mutable struct two_phase_modelType{T1, T2, M}
    m1::Type{T1}
    m2::Type{T2}
    mix::M
end

"""
    two_phase_model(ϕ, m1, m2, mix)

Rock physics model to combine two phases, usually constructed through `two_phase_modelType`[@ref]

## Arguments
 - `ϕ` : Vol. fraction of the **second** phase
 - `m1` : model corresponding to phase 1
 - `m2` : model corresponding to phase 2
 - `mix` : mixing type, available options are `HS_1962_plus()`, `HS1962_minus`, `MAL(m)`

## Usage

```julia
m = two_phase_modelType(SEO3, Ni2011, HS1962_plus())
ps_nt = ps_nt = (; T=[800.0f0, 1000.0f0] .+ 273, P=3.0f0, ρ=3300.0f0, Ch2o_m=1000.0f0, ϕ=0.1f0)
model = m(ps_nt)

resp = forward(model)
```
"""
mutable struct two_phase_model{V, T1, T2, M} # <: AbstractRockphyModel
    ϕ::V
    m1::T1
    m2::T2
    mix::M
end

two_phase_modelType(m1) = m1
two_phase_modelType(m1, m::phase_mixing) = m1

function (model::two_phase_modelType)(ps::NamedTuple)
    mix = model.mix
    ϕ = ps.ϕ
    v1 = from_nt(getproperty(model, :m1), ps)
    v2 = from_nt(getproperty(model, :m2), ps)
    return two_phase_model(ϕ, v1, v2, mix)
end

function forward(model::two_phase_model{V, T1, T2, M},
        p) where {V, M, T1 <: AbstractCondModel, T2 <: AbstractCondModel}
    σ1 = MT.forward(model.m1, []).σ
    σ2 = MT.forward(model.m2, []).σ

    @. σ1 = exp10(σ1)
    @. σ2 = exp10(σ2)

    σ = broadcast(
        (sig1, sig2, phi) -> MT.mix_models([sig1, sig2], phi, model.mix), σ1, σ2, model.ϕ)

    return RockphyCond(log10.(σ))
end

function forward(model::two_phase_model{V, T1, T2, M}, p,
        params) where {V, M, T1 <: AbstractCondModel, T2 <: AbstractCondModel}
    σ1 = MT.forward(model.m1, [], params.m1).σ
    σ2 = MT.forward(model.m2, [], params.m2).σ

    @. σ1 = exp10(σ1)
    @. σ2 = exp10(σ2)

    σ = broadcast(
        (sig1, sig2, phi) -> MT.mix_models([sig1, sig2], phi, model.mix), σ1, σ2, model.ϕ)

    return RockphyCond(log10.(σ))
end

function default_params(::Type{two_phase_model{V, T1, T2, M}}) where {V, T1, T2, M}
    (; zip([:m1, :m2], [default_params(T1), default_params(T2)])...)
end

# following is needed for combine_models

function from_nt(m::Type{T}, nt::NamedTuple) where {T <: two_phase_modelType}

    ϕ = nt.ϕ
    m1 = m.types[1].parameters[1]
    m2 = m.types[2].parameters[1]
    mix = m.types[3]

    model1 = MT.from_nt(m1, nt)
    model2 = MT.from_nt(m2, nt)

    return two_phase_model(ϕ, model1, model2, mix())
end

function from_nt(m::Type{T}, nt::NamedTuple) where {T <: two_phase_model}

    ϕ = nt.ϕ
    m1 = T.parameters[2]
    m2 = T.parameters[3]
    mix = T.parameters[4]

    model1 = MT.from_nt(m1, nt)
    model2 = MT.from_nt(m2, nt)

    return two_phase_model(ϕ, model1, model2, mix())
end

#= ==============================================================================
multi-phase (stochastic inverse looks hard with this) : Experimental
=#

mutable struct construct_model_multi_phase2{T1, T2, T3, T4, T5, M}
    m1::Type{T1}
    m2::Type{T2}
    m3::Type{T3}
    m4::Type{T4}
    m5::Type{T5}
    mix::M
end

construct_model_multi_phase2(m1) = m1
construct_model_multi_phase2(m1, m::phase_mixing) = m1

function construct_model_multi_phase2(m1, m2, m::phase_mixing)
    construct_model_multi_phase2(m1, m2, Nothing, Nothing, Nothing, m)
end
function construct_model_multi_phase2(m1, m2, m3, m::phase_mixing)
    construct_model_multi_phase2(m1, m2, m3, Nothing, Nothing, m)
end
function construct_model_multi_phase2(m1, m2, m3, m4, m::phase_mixing)
    construct_model_multi_phase2(m1, m2, m3, m4, Nothing, m)
end

# @inferred construct_model_multi_phase2(SEO3, Ni2011, HS1962_plus())
# m1 = construct_model_multi_phase2(SEO3, Ni2011, HS1962_plus())

# function from_nt(m::Type{T}, nt::NamedTuple) where T<:construct_model_multi_phase

#     # fnames = fieldnames(T)
#     ϕ = getproperty(ps_nt, :ϕ)
#     m1 = m.types[1].parameters[1]
#     m2 = m.types[2].parameters[1]
#     m3 = m.types[3].parameters[1]
#     m4 = m.types[4].parameters[1]
#     m5 = m.types[5].parameters[1]
#     mix = m.types[4]

#     model1 = MT.from_nt(m1, ps_nt)
#     model2 = MT.from_nt(m2, ps_nt)

#     return model_multi_phase(ϕ, model1, model2, mix())

# end

mutable struct model_multi_phase2{V, T1, T2, T3, T4, T5, M}
    ϕ::V
    m1::T1
    m2::T2
    m3::T3
    m4::T4
    m5::T5
    mix::M
end

function rearrange_ϕ(ϕ, model::construct_model_multi_phase2)
    @assert sum(ϕ)≤1 "Σϕᵢ = $(sum(ϕ)) should be ≤ 1."

    fnames = propertynames(model)[1:(end - 1)]
    fnames = filter(f -> getfield(model, f) !== Nothing, fnames)
    msg = """
    Vol. frac of last component is defined automatically once the others are defined. 
    Make sure the length of porosity parameter `ϕ` is one less than the length of number of components.
    length of ϕ = $(length(ϕ))
    length of components = $(length(fnames))
    Check out the relevant documentation.
    """

    @assert length(ϕ)==(length(fnames) - 1) msg

    c = length(ϕ)
    ϕ_vec = zeros(eltype(ϕ), 5)

    ϕ_vec[1:c] .= ϕ
    ϕ_vec[c + 1] = 1 - sum(ϕ)
    return ϕ_vec
end

function (model::construct_model_multi_phase2)(ps::NamedTuple)
    pnames = propertynames(model)
    mix = getfield(model, pnames[end])
    # ϕ = getfield(ps, :ϕ)

    # make ϕ according to different number of phases
    ϕ_vec = rearrange_ϕ(ps.ϕ, model)

    # pnames = pnames[2:3]
    v1 = from_nt(getproperty(model, pnames[1]), ps)
    v2 = from_nt(getproperty(model, pnames[2]), ps)
    v3 = from_nt(getproperty(model, pnames[3]), ps)
    v4 = from_nt(getproperty(model, pnames[4]), ps)
    v5 = from_nt(getproperty(model, pnames[5]), ps)
    return model_multi_phase2(ϕ_vec, v1, v2, v3, v4, v5, mix)
end

# model = m1(ps_nt)

function MT.forward(model::model_multi_phase2{V, T1, T2, T3, T4, T5, M},
        p) where {V, M <: Union{HS1962_minus, HS1962_plus, MAL},
        T1 <: AbstractCondModel, T2, T3, T4, T5}
    σ1 = MT.forward(model.m1, []).σ
    σ2 = MT.forward(model.m2, []).σ

    @. σ1 = exp10(σ1)
    @. σ2 = exp10(σ2)

    σ = broadcast(
        (sig1, sig2, phi) -> MT.mix_models([sig1, sig2], phi, model.mix), σ1, σ2, model.ϕ)

    return RockphyCond(log10.(σ))
end
