mutable struct construct_model_2phase{T1, T2, M}
    m1::Type{T1}
    m2::Type{T2}
    mix::M
end

mutable struct model_2phase{V, T1, T2, M}
    ϕ::V
    m1::T1
    m2::T2
    mix::M
end

construct_model_2phase(m1) = m1
construct_model_2phase(m1, m::phase_mixing) = m1
# construct_model_2phase(m1, m2, m::phase_mixing) = construct_model_2phase(m1, m2, m)

function (model::construct_model_2phase)(ps::NamedTuple)
    mix = model.mix
    ϕ = ps.ϕ
    v1 = from_nt(getproperty(model, :m1), ps)
    v2 = from_nt(getproperty(model, :m2), ps)
    return model_2phase(ϕ, v1, v2, mix)
end

function forward(model::model_2phase{V, T1, T2, M},
        p; params = default_params(Val{model_2phase{V, T1, T2, M}}())) where {V, M, T1 <: AbstractCondModel, T2 <: AbstractCondModel}
    σ1 = MT.forward(model.m1, [], params = params.m1).σ
    σ2 = MT.forward(model.m2, [], params = params.m2).σ

    @. σ1 = exp10(σ1)
    @. σ2 = exp10(σ2)

    σ = broadcast(
        (sig1, sig2, phi) -> MT.mix_models([sig1, sig2], phi, model.mix), σ1, σ2, model.ϕ)

    return RockphyCond(log10.(σ))
end

# following is needed for combine_models
function from_nt(m::Type{T}, ps_nt::NamedTuple) where {T <: construct_model_2phase}

    # fnames = fieldnames(T)
    ϕ = getproperty(ps_nt, :ϕ)
    m1 = m.types[1].parameters[1]
    m2 = m.types[2].parameters[1]
    mix = m.types[3]

    model1 = MT.from_nt(m1, ps_nt)
    model2 = MT.from_nt(m2, ps_nt)

    return model_2phase(ϕ, model1, model2, mix())
end

function default_params(::Val{model_2phase{V, T1, T2, M}}) where {V, T1, T2, M}
    (; zip(
        [:m1, :m2],
        [default_params(Val{T1}()), default_params(Val{T2}())]
    )...)
end

#= ==============================================================================
multi-phase (stochastic inverse would be hard with this)
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
