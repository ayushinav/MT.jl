mutable struct construct_model_multi_rp{T1, T2, T3, T4}
    cond::Type{T1}
    elastic::Type{T2}
    visc::Type{T3}
    anelastic::Type{T4}
end

mutable struct model_multi_rp{T1, T2, T3, T4}
    cond::T1
    elastic::T2
    visc::T3
    anelastic::T4
end

function (model::construct_model_multi_rp)(ps::NamedTuple)
    v = map((x) -> from_nt(getproperty(model, x), ps), propertynames(model))
    return model_multi_rp(v...)
end

@generated function forward(model::model_multi_rp{T1, T2, T3, T4}, p) where {T1, T2, T3, T4}
    fnames = fieldnames(model_multi_rp)
    args = [:(to_nt(MT.forward(getfield(model, $(QuoteNode(f))), []))) for f in fnames]
    return :(merge($(args...)))
end
