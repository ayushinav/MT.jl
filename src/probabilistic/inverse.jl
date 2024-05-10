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
* `r_obs`: `response` that needs to inverted for
* `err_resp`: `response` variable containing the errors associated with observed response
* `vars`: variables that need to be passed into the `forward` function along with `model` to generate a `response`
* `alg_cache`: to tell the compiler what type of stochastic inversion method is to be used
* `trans_utils`: A named tuple containing `transform_utils` for the fields of model that need to be scaled/modified. If not provided for any `model` field, the field won't be modified.
"""
function stochastic_inverse(
    r_obs::response,
    err_resp::response,
    vars,
    alg_cache::mcmc_cache;
    trans_utils::NamedTuple = (m = log_tf, h = lin_tf)
    )

    model_fields = [];
    modelD = [];
    const_data = []; #model([], []);

    # segregate the constants and the Distribution parts of the alg_cache

    for k in fieldnames(typeof(alg_cache.apriori))
        if typeof(getfield(alg_cache.apriori, k)) <: Distribution
            push!(model_fields, k);
            push!(modelD, getfield(alg_cache.apriori, k));
            push!(const_data, rand(getfield(alg_cache.apriori, k)))
        else
            push!(const_data, getfield(alg_cache.apriori, k))
        end
    end

    response_fields = Symbol.([])
    for k in fieldnames(typeof(alg_cache.likelihood))
        if typeof(getfield(alg_cache.likelihood, k)) <: Function
            push!(response_fields, k)
        end
    end
    
    # putting trans_utils together for all the fields

    trans_utils_arr = [];
    for k ∈ fieldnames(typeof(alg_cache.apriori))
        if k in keys(trans_utils)
            push!(trans_utils_arr, trans_utils[k])
        else
            push!(trans_utils_arr, lin_tf)
        end
    end

    transf_utils = (;zip(
        [fieldnames(typeof(alg_cache.apriori))...], trans_utils_arr)...);

    # make model and modelDistribution

    m_sample = model(const_data...);
    m₀ = (;zip(
        [fieldnames(typeof(m_sample))...], 
        [getfield(m_sample, k) for k ∈ fieldnames(typeof(m_sample))]
    )...);

    robs = (;zip(
        [fieldnames(typeof(r_obs))...], 
        [getfield(r_obs, k) for k ∈ fieldnames(typeof(r_obs))]
    )...);
    
    mcmc_model = mcmc_turing(
        # m₀, # ::NamedTuple
        m_sample, # ::model
        vars,
        robs, # ::NamedTuple
        err_resp, # ::response
        alg_cache.apriori, # ::NamedTuple
        alg_cache.likelihood, # ::responseDistribution
        response_fields = Symbol.(response_fields),
        model_fields = Symbol.(model_fields),
        trans_utils = transf_utils
    )

    chains = sample(mcmc_model, alg_cache.sampler, alg_cache.n_samples)
    return chains;
end
