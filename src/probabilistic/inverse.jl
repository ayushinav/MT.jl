"""
    stochastic_inverse(
        r_obs::response,
        err_resp::response,
        vars,
        alg_cache::mcmc_cache;
        model_trans_utils::NamedTuple = (m = lin_tf, h = lin_tf)
        )

function to perform sampling

### Returns

    `AbstractMCMC.jl Chain` containing the vectors used for performing samping. Note that all the variables will be named `var_inf`, for variable inferred. 
    The vector in each sample will be all the fields of the model being inferred on concatenated together.

### Variables

  - `r_obs`: `response` that needs to inverted for
  - `err_resp`: `response` variable containing the errors associated with observed response
  - `vars`: variables that need to be passed into the `forward` function along with `model` to generate a `response`
  - `alg_cache`: to tell the compiler what type of stochastic inversion method is to be used
  - `model_trans_utils`: A named tuple containing `transform_utils` for the fields of model that need to be scaled/modified. If not provided for any `model` field, the field won't be modified.
"""
function stochastic_inverse(r_obs::resp1, err_resp::resp2, vars, alg_cache::mcmc_cache;
        model_trans_utils::NamedTuple=(m=lin_tf, h=lin_tf), # need to take care of this
        response_trans_utils::NamedTuple=(; ρₐ=lin_tf, ϕ=lin_tf), params=(;),
        kwargs...) where {resp1 <: AbstractResponse, resp2 <: AbstractResponse}
    model_fields = Symbol[]
    # modelD = []
    const_data = []

    # segregate the constants and the Distribution parts of the alg_cache

    apriori = to_dist_nt(alg_cache.apriori)

    for k in keys(apriori)
        if typeof(getfield(apriori, k)) <: Distribution
            push!(model_fields, k)
            push!(const_data, rand(getfield(apriori, k)))
        else
            push!(const_data, getfield(apriori, k))
        end
    end

    response_fields = Symbol[]
    likelihood = to_dist_nt(alg_cache.likelihood)
    for k in keys(likelihood) # similarly, here it will be propertynames for likelihood being a NamedTuple
        if typeof(getfield(likelihood, k)) <: Function
            push!(response_fields, k)
        end
    end

    # putting trans_utils together for all the fields

    trans_utils_arr = []
    for k in keys(apriori)
        if k in keys(model_trans_utils)
            push!(trans_utils_arr, model_trans_utils[k])
        else
            push!(trans_utils_arr, lin_tf)
        end
    end

    transf_utils = (; zip(keys(apriori), trans_utils_arr)...) # NamedTuple for trans_utils and defaults

    m_type = sample_type(alg_cache.apriori)

    if isempty(params)
        params = default_params(m_type)
    end

    # @show to_resp_nt(r_obs)
    # @show to_resp_nt(err_resp)
    msg = """
    variables to be inferred : $(model_fields)
    variables used for inference : $(response_fields)
    """
    @info msg

    mcmc_model = mcmc_turing(m_type, const_data, vars, to_resp_nt(r_obs), # ::NamedTuple
        to_resp_nt(err_resp), # ::response
        apriori, # ::NamedTuple
        likelihood, # ::responseDistribution
        params; response_fields=Symbol.(response_fields),
        model_fields=Symbol.(model_fields), model_trans_utils=transf_utils,
        response_trans_utils=response_trans_utils)

    if typeof(alg_cache.sampler) <: Turing.AdvancedVI.VariationalInference
        return vi(mcmc_model, alg_cache.sampler)
    else
        return Turing.sample(
            mcmc_model, alg_cache.sampler, alg_cache.n_samples; verbose=false, kwargs...)
    end
end
