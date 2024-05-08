
function inverse!( 
    r_obs::response,
    vars,
    alg_cache::mcmc_cache;
    trans_utils::NamedTuple = (m = log_tf, h = lin_tf), # change this to a NamedTuple
    verbose = false)

    # test if model size is compatible with model distribution, same with response and response distribution

    model_fields = [];
    modelD = [];
    const_data = [];

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
        if typeof(getfield(alg_cache.likelihood, k)) <: Distribution
            push!(response_fields, k)
        end
    end


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
    
    mDist = (;zip(model_fields, modelD)...);

    mcmc_model = mcmc_turing(
        m₀, # ::NamedTuple
        m_sample, # ::model
        vars,
        robs, # ::NamedTuple
        mDist, # ::NamedTuple
        alg_cache.likelihood, # ::responseDistribution2;
        response_fields = Symbol.(response_fields), #::Vector{Symbol}= [k for k ∈ fieldnames(typeof(rDist))],
        model_fields = Symbol.(model_fields),
        trans_utils = transf_utils
    )

    chains = sample(mcmc_model, alg_cache.sampler, alg_cache.n_samples)
    return chains;
end