mutable struct multi_rp_response{T1, T2, T3, T4} <: AbstractRockphyResponse
    cond::T1
    elastic::T2
    visc::T3
    anelastic::T4
end

mutable struct construct_model_multi_rp{T1, T2, T3, T4}
    cond::Type{T1}
    elastic::Type{T2}
    visc::Type{T3}
    anelastic::Type{T4}
end

mutable struct model_multi_rp{T1, T2, T3, T4} # <: AbstractRockphyModel
    cond::T1
    elastic::T2
    visc::T3
    anelastic::T4
end

function (model::construct_model_multi_rp)(ps::NamedTuple)
    v = map((x) -> from_nt(getproperty(model, x), ps), propertynames(model))
    return model_multi_rp(v...)
end

function forward(model::model_multi_rp{T1, T2, T3, T4}, p, params) where {T1, T2, T3, T4}
    resp_cond = forward(model.cond, [], params.cond)
    resp_elastic = forward(model.elastic, [], params.elastic)
    resp_visc = forward(model.visc, [], params.visc)
    resp_anelastic = forward(model.anelastic, [], params.anelastic)

    # (; to_nt(resp_cond)..., to_nt(resp_elastic)..., to_nt(resp_visc)..., to_nt(resp_anelastic)...)
    multi_rp_response(resp_cond, resp_elastic, resp_visc, resp_anelastic)
end

function forward(model::model_multi_rp{T1, T2, T3, T4}, p) where {T1, T2, T3, T4}
    resp_cond = forward(model.cond, [])
    resp_elastic = forward(model.elastic, [])
    resp_visc = forward(model.visc, [])
    resp_anelastic = forward(model.anelastic, [])

    # (; to_nt(resp_cond)..., to_nt(resp_elastic)..., to_nt(resp_visc)..., to_nt(resp_anelastic)...)
    multi_rp_response(resp_cond, resp_elastic, resp_visc, resp_anelastic)
end

function default_params(::Type{model_multi_rp{T1, T2, T3, T4}}) where {T1, T2, T3, T4}
    (;
        zip([:cond, :elastic, :visc, :anelastic],
            [default_params(T1), default_params(T2),
                default_params(T3), default_params(T4)])...)
end
