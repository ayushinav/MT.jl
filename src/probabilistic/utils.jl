mutable struct mcmc_cache
    apriori::modelDistribution
    likelihood::responseDistribution
    n_samples::Int
    sampler
end

@model function mcmc_turing(
    m₀::NamedTuple, # let's see
    m_sample::model,
    vars,
    r_obs::NamedTuple,
    mDist::NamedTuple,
    rDist::responseDistribution;
    response_fields::Vector{Symbol}= [k for k ∈ fieldnames(typeof(rDist))],
    model_fields::Vector{Symbol}= [k for k ∈ fieldnames(typeof(mDist))],
    trans_utils::NamedTuple = (m = log_tf, h = lin_tf)
    ) 
    
    for k in model_fields
        kk ~ mDist[k]
        broadcast!(getproperty(trans_utils[k], :itf), m₀[k], kk)
    end

    r_sample = forward(m_sample, vars);

    # formulate rDist using r_sample
    update_responseDistribution!(rDist, r_sample)

    for k in response_fields
        r_obs[k] ~ getfield(rDist, k)
    end
end