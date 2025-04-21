mutable struct construct_model_2_phase{V, T1, T2, M}
    ϕ::V
    m1::Type{T1}
    m2::Type{T2}
    mix::M
end

mutable struct model_2_phase{V, T1, T2, M}
    ϕ::V
    m1::T1
    m2::T2
    mix::M
end

construct_model_2_phase(ϕ::V, m1) where V <: AbstractVector = m1
construct_model_2_phase(ϕ::V, m1, m::phase_mixing) where V <: AbstractVector = m1
construct_model_2_phase(ϕ::V, m1, m2, m::phase_mixing) where V <: AbstractVector = construct_model_2_phase(ϕ, m1, m2, Nothing, Nothing, Nothing, m)

# @inferred construct_model_2_phase4(0.1f0, SEO3, Ni2011, HS1962_plus())
# m1 = construct_model_2_phase4(0.1f0, SEO3, Ni2011, HS1962_plus())


function (model::construct_model_2_phase)(ps::NamedTuple)

    pnames = propertynames(model)
    mix = getfield(model, pnames[end])
    ϕ = getfield(model, pnames[1])
    # pnames = pnames[2:3]
    v1 = from_nt(getproperty(model, pnames[2]), ps)
    v2 = from_nt(getproperty(model, pnames[3]), ps)
    return model_2_phase4(ϕ, v1, v2, mix)

end

# ps_nt = (; T = [1200f0, 1400f0], P = 3f0, ρ = 3300f0, Ch2o_m = 100f0)

# @inferred m1(ps_nt)
# model = m1(ps_nt)

function forward(model::model_2_phase4{V, T1, T2, M}) where {V, M, T1 <: AbstractCondModel, T2 <: AbstractCondModel}
    σ1 = MT.forward(model.m1, []).σ
    σ2 = MT.forward(model.m2, []).σ

    @. σ1 = exp10(σ1)
    @. σ2 = exp10(σ2)

    σ = broadcast((sig1, sig2, phi) -> MT.mix_models([sig1, sig2], phi, model.mix), σ1, σ2, model.ϕ)

    return RockphyCond(log10.(σ))

end

# @inferred forward(model)
forward(model)


