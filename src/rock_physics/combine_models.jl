
mutable struct model_multi_rp{T2}
    # rp_type::T1
    ps::T2 # vector of parameters 
    p_names::Vector{Symbol} # Vector of symbols telling the parameters in vector
    model_list::Vector{<:Any}
    params
    resp_fields
    # Ch2o::T4
    # h2o_part::T5
end

mutable struct construct_model_multi_rp
    model_list
    resp_fields

    function construct_model_multi_rp(model_list, resp_fields)
        resp_list = forward.(model_list)

        var_list = vcat([[fieldnames(im)...] for im in resp_list]...)

        f_ = reduce(&, [ir ∈ var_list for ir in resp_fields])

        @assert f_ === true """the responses asked for cannot be obtained form the models provided. \n
            parameters asked for : $resp_fields \n
            all parameters models can output : $var_list \n 
        """

        return new(model_list, resp_fields)
    end
end

function construct_model_multi_rp(model_list)
    resp_list = forward.(model_list)
    var_list = unique(vcat([[fieldnames(im)...] for im in resp_list]...))

    return construct_model_multi_rp(model_list, var_list)
end

function (m::construct_model_multi_rp)(;params = (;), p...)
    var_list = vcat([[fieldnames(ir)...] for ir in m.model_list]...)
    unique!(var_list)
    
    p_names = [keys(p)...]
    ps = [p[k] for k in p_names]
    
    params_vec = [(Symbol(im) ∈ keys(params)) ? (getfield(params, Symbol(im))) : MT.default_params(Val{im}()) for im in m.model_list]
    params_nt =  (;zip(Symbol.(m.model_list), (params_vec))...)

    f_ = reduce(&, [ir ∈ p_names for ir in var_list])
    @assert f_ == true """all the variables required by models are not included in `p_names`. \n
      parameters mentioned : $p_names \n
      all parameters required : $var_list \n 
      """ # can have a better message here (p_names belonging to corresponding rp_types)

    resp_list = [forward(im) for im in m.model_list]
    var_list = vcat([[fieldnames(im)...] for im in resp_list]...)

    return model_multi_rp(ps, p_names, m.model_list, params_nt, m.resp_fields)

end

rp_property_names(x) = fieldnames(x)
rp_property_names(x::T) where {T <: MT.model_multiphase} = x.p_names

function forward(m::model, p) where {model <: model_multi_rp}
    resps = []; #zeros(eltype(m.ps), length(m.model_list))

    ps = (; zip(m.p_names, m.ps)...)

    # water partition code
        
    for i in eachindex(m.model_list)
        var_list = [(rp_property_names(m.model_list[i]))...]
        @show m.model_list[i](ps[var_list]...)
        resp = forward(m.model_list[i](ps[var_list]...), []; params = getfield(m.params, Symbol(m.model_list[i])))
        # resps[i] = getfield(resp, first(m.resp_fields))
        push!(resps, resp)
    end
    
    var_list = vcat([[propertynames(ir)...] for ir in resps]...)
    vals_list = vcat([[[getfield(ir, im) for im in propertynames(ir)]...] for ir in resps]...)
    
    resp_net = (; zip(var_list, vals_list)...) # or maybe output a vector of symbols and a vector of values

    return resp_net
end
