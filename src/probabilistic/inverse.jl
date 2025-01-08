"""
    stochastic_inverse(
        r_obs::response,
        err_resp::response,
        vars,
        alg_cache::mcmc_cache;
        trans_utils::NamedTuple = (m = log_tf, h = lin_tf)
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
  - `trans_utils`: A named tuple containing `transform_utils` for the fields of model that need to be scaled/modified. If not provided for any `model` field, the field won't be modified.
"""
function stochastic_inverse(r_obs::resp1, err_resp::resp2, vars, alg_cache::mcmc_cache;
        trans_utils::NamedTuple=(m=log_tf, h=lin_tf), # need to take care of this
        kwargs...) where {resp1 <: AbstractResponse, resp2 <: AbstractResponse}
    model_fields = []
    modelD = []
    const_data = [] #model([], []);

    # segregate the constants and the Distribution parts of the alg_cache

    for k in fieldnames(typeof(alg_cache.apriori)) # make it properynames and make alg_cache.apriori a NamedTuple
        if typeof(getfield(alg_cache.apriori, k)) <: Distribution # getfield will be replaced by getproperty
            push!(model_fields, k)
            push!(modelD, getfield(alg_cache.apriori, k))
            push!(const_data, rand(getfield(alg_cache.apriori, k)))
        else
            push!(const_data, getfield(alg_cache.apriori, k))
        end
    end

    response_fields = Symbol.([])
    for k in fieldnames(typeof(alg_cache.likelihood)) # similarly, here it will be propertynames for likelihood being a NamedTuple
        if typeof(getfield(alg_cache.likelihood, k)) <: Function
            push!(response_fields, k)
            # @show k
        end
        # @show k, typeof(k)
    end

    # putting trans_utils together for all the fields

    trans_utils_arr = []
    for k in fieldnames(typeof(alg_cache.apriori))
        if k in keys(trans_utils)
            push!(trans_utils_arr, trans_utils[k])
        else
            push!(trans_utils_arr, lin_tf)
        end
    end

    transf_utils = (; zip([fieldnames(typeof(alg_cache.apriori))...], trans_utils_arr)...) # NamedTuple for trans_utils and defaults

    # make model and modelDistribution
    # @show const_data
    # @show MT.inverse(r_obs, abstract = true)
    # @show response_fields

    m_sample = MT.inverse(r_obs; abstract=true)(const_data...)

    robs = (;
        zip([fieldnames(typeof(r_obs))...],
            [getfield(r_obs, k) for k in fieldnames(typeof(r_obs))])...)

    mcmc_model = mcmc_turing(m_sample, const_data, vars, robs, # ::NamedTuple
        err_resp, # ::response
        alg_cache.apriori, # ::NamedTuple
        alg_cache.likelihood; # ::responseDistribution
        response_fields=Symbol.(response_fields),
        model_fields=Symbol.(model_fields), trans_utils=transf_utils)

    # print("MODEL WORKS?")
    # @show typeof(sampler) <: Turing.AdvancedVI.VariationalInference
    @show kwargs

    if typeof(alg_cache.sampler) <: Turing.AdvancedVI.VariationalInference
        return vi(mcmc_model, alg_cache.sampler)
    else
        return Turing.sample(
            mcmc_model, alg_cache.sampler, alg_cache.n_samples; verbose=false, kwargs...)
    end
end

import Turing:AbstractChains

