mutable struct multi_rp_responseDistribution{T1, T2, T3, T4} <:
               AbstractRockphyResponseDistribution
    cond::T1
    elastic::T2
    visc::T3
    anelastic::T4
end

mutable struct construct_model_multi_rpDistribution{T1, T2, T3, T4}
    cond::Type{T1}
    elastic::Type{T2}
    visc::Type{T3}
    anelastic::Type{T4}
end

mutable struct model_multi_rpDistribution{T1, T2, T3, T4} <:
               AbstractRockphyModelDistribution
    cond::T1
    elastic::T2
    visc::T3
    anelastic::T4
end

function (model::construct_model_multi_rpDistribution)(ps::NamedTuple)
    v = map((x) -> from_nt(getproperty(model, x), ps), propertynames(model))
    return model_multi_rpDistribution(v...)
end
