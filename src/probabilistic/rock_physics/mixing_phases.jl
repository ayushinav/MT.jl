mutable struct construct_model_2phaseDistribution{T1, T2, M}
    m1::Type{T1}
    m2::Type{T2}
    mix::M
end

mutable struct model_2phaseDistribution{V, T1, T2, M} <: AbstractRockphyModelDistribution
    ϕ::V
    m1::T1
    m2::T2
    mix::M
end

construct_model_2phaseDistribution(m1) = m1
construct_model_2phaseDistribution(m1, m::phase_mixing) = m1

function (model::construct_model_2phaseDistribution)(ps::NamedTuple)
    mix = model.mix
    ϕ = ps.ϕ

    v1 = from_nt(model.m1, ps)

    v2 = from_nt(model.m2, ps)

    return model_2phaseDistribution(ϕ, v1, v2, mix)
end
