mutable struct mcmc_cache
    apriori::modelDistribution
    likelihood::responseDistribution
    n_samples::Int
    sampler
end

@model function mcmc_turing(
    m_sample::model,
    vars,
    r_obs::NamedTuple,
    err_resp::response,
    mDist::modelDistribution,
    rDist::responseDistribution;
    response_fields::Vector{Symbol}= [k for k ∈ fieldnames(typeof(rDist))],
    model_fields::Vector{Symbol}= [k for k ∈ fieldnames(typeof(mDist))],
    trans_utils::NamedTuple = (m = log_tf, h = lin_tf)
    )

    m₀ = [];
    for k ∈ propertynames(mDist)
        if k in model_fields
            kk ~ getproperty(mDist, k)
            push!(m₀, broadcast(getproperty(trans_utils[k], :itf), kk));
        else
            push!(m₀, getfield(m_sample, k));   
        end
    end

    m_sample = model(m₀...);
    r_sample = forward(m_sample, vars);

    for k in response_fields
        r_obs[k] ~ getfield(rDist, k)(getfield(r_sample, k), getfield(err_resp, k))
    end
end