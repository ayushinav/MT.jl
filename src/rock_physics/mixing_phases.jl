
mutable struct model_multiphase{T2, T3} <: MT.AbstractRockphyModel
    ps::T2 # vector of parameters 
    p_names::Vector{Symbol} # Vector of symbols telling the parameters in vector
    ϕ::T3 # phase ratios
    model_list::Vector{<:Any}
    mixing_type
    params
    resp_fields
    # Ch2o::T4
    # h2o_part::T5
end

mutable struct consturct_model_multiphase
    model_list
    mixing_type
    resp_fields
    # p_names
end

function consturct_model_multiphase(model_list, mixing_type)
    if !(reduce(&, supertype.(model_list) .<: supertype(first(model_list))))
        @error "all the model names should be of the same supertype. If they are of different types, provide a `resp_fields` argument"
    end

    rp_type = forward(first(model_list))
    resp_fields = [fieldnames(rp_type)...]

    if typeof(mixing_type) <: Union{HS1962_plus, HS1962_minus, MAL}
        @assert length(model_list)==2 "`$mixing_type` model allows for only 2 models to mix, the first one being the solid and the second melt"

    elseif typeof(mixing_type) <: single_phase
        @assert length(model_list)==1 "single phase models require only one model in `model_list`"
    end

    return consturct_model_multiphase(model_list, mixing_type, resp_fields)
end

model_mixed = consturct_model_multiphase([SEO3, Ni2011], HS1962_plus())

function (m::consturct_model_multiphase)(; ϕ=0.0f0, params=(;), p...)
    var_list = vcat([[fieldnames(ir)...] for ir in m.model_list]...)
    unique!(var_list)
    p_names = [keys(p)...]
    ps = [p[k] for k in p_names]

    params_vec = [(Symbol(im) ∈ keys(params)) ? (getfield(params, Symbol(im))) :
                  MT.default_params(Val{im}()) for im in m.model_list]
    params_nt = (; zip(Symbol.(m.model_list), (params_vec))...)

    f_ = reduce(&, [ir ∈ p_names for ir in var_list])
    @assert f_==true """all the variables required by models are not included in `p_names`. \n
    parameters mentioned : $p_names \n
    all parameters required : $var_list \n 
    """ # can have a better message here (p_names belonging to corresponding rp_types)

    if typeof(m.mixing_type) <: Union{HS1962_plus, HS1962_minus, MAL}
        @assert length(ϕ)==1 "only the fraction of the second component, ie melt, $(model_list[2]) is needed"

        return model_multiphase(
            ps, p_names, ϕ, m.model_list, [m.mixing_type], params_nt, m.resp_fields)

    elseif typeof(m.mixing_type) <: single_phase
        @assert length(ϕ) == 1&&first(ϕ) == 1 "single phase models will have total fraction for the single phase"
        return model_multiphase(
            ps, p_names, ϕ, m.model_list, [m.mixing_type], params_nt, m.resp_fields)

    else
        @assert length(model_list)==(length(ϕ) + 1) """each phase and it's vol fraction should be provided, except for the last one, 
        where the vol fraction for the last one is obtained from the remaining variables 1 - ∑ϕ"""
        return model_multiphase(
            ps, p_names, ϕ, m.model_list, [m.mixing_type], params_nt, m.resp_fields)
    end
end

function forward(m::model, p) where {model <: mixers}
    resps = [] #zeros(eltype(m.ps), length(m.model_list))
    # resps_ = []

    ps = (; zip(m.p_names, m.ps)...)

    # water partition code

    for i in eachindex(m.model_list)
        var_list = [(fieldnames(m.model_list[i]))...]
        resp = MT.forward(m.model_list[i](ps[var_list]...), [];
            params=getfield(m.params, Symbol(m.model_list[i])))
        # resps[i] = getfield(resp, first(m.resp_fields))
        push!(resps, resp)
    end

    var_list = [(fieldnames(first(m.model_list)))...]
    resp_net = MT.forward(first(m.model_list)(ps[var_list]...), [];
        params=getfield(m.params, Symbol(first(m.model_list))))

    for k in m.resp_fields
        r_net = MT.mix_models(getfield.(resps, k), m.ϕ, first(m.mixing_type))
        setfield!(resp_net, k, r_net)
    end

    return resp_net
end

function mix_models(σs, ϕ, ::HS1962_plus)
    σ_max = 10.0f0^maximum(σs)
    σ_min = 10.0f0^minimum(σs)
    phi = first(ϕ)

    num = 3 * (1 - phi) * (σ_max - σ_min) # numerator
    den = 3 * σ_max - phi * (σ_max - σ_min) # denominator
    esig = σ_max * (1 - (num / den))

    return log10(esig)
end

function mix_models(σs, ϕ, ::HS1962_minus)
    σ_max = 10.0f0^maximum(σs)
    σ_min = 10.0f0^minimum(σs)
    phi = first(ϕ)

    num = 3 * (phi) * (σ_max - σ_min) # numerator
    den = 3 * σ_min + (1 - phi) * (σ_max - σ_min) # denominator
    esig = σ_min * (1 + (num / den))

    return log10(esig)
end

function mix_models(σs, ϕ, mal::MAL)
    σ_fluid = 10.0f0^(σs[2])
    σ_matrix = 10.0f0^(σs[1])

    phi = first(ϕ)
    sig = σ_fluid

    if phi < 1
        p = log10(1 - phi^mal.m) * inv(log10(1 - phi))
        sig = σ_fluid * phi^mal.m + σ_matrix * (1 - phi)^p
    end

    return log10(sig)
end

function mix_models(σs, ϕ, ::single_phase)
    @assert length(σs) == 1
    return first(σs)
end
