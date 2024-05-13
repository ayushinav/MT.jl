"""
    mutable struct mcmc_cache
        apriori::modelDistribution
        likelihood::responseDistribution
        n_samples::Int
        sampler
    end
placeholder to store 
* apriori in the form of [`modelDistribution`](@ref modelDistribution)
* likelihood in the form of [`responseDistribution`](@ref responseDistribution)
* number of samples to obtain in `n_samples`
* `Turing.jl` sampler to be used in `sampler`
"""
mutable struct mcmc_cache
    apriori::MTModelDistribution
    likelihood::MTResponseDistribution
    n_samples::Int
    sampler
end


"""
    @model function mcmc_turing(
        vars,
        r_obs::NamedTuple,
        err_resp::response,
        mDist::modelDistribution,
        rDist::responseDistribution;
        response_fields::Vector{Symbol}= [k for k ∈ fieldnames(typeof(rDist))],
        model_fields::Vector{Symbol}= [k for k ∈ fieldnames(typeof(mDist))],
        trans_utils::NamedTuple = (m = log_tf, h = lin_tf)
        )
makes a `Turing.jl` model to perform MCMC sampling

### Variables:
* `vars`: variables that need to be passed into the `forward` function along with `model` to generate a `response`
* `r_obs`: named tuple containing the observed data, with the same keys as the fields in the corresponding `response`
* `err_resp`: `response` variable that contains the errors
* `mDist`: [`modelDistribution`](@ref modelDistribution) contains the apriori information
* `rDist`: [`responseDistribution`](@ref responseDistribution) contains the likelihood information

### Keyword/optional arguments
* `response_fields`:  which fields in `response` to invert for
* `model_fields`: fields in `model` to draw inference on
* `trans_utils`: to transform the model field variables to and from computational (inference) domain
"""
@model function mcmc_turing(
    m_sample::model,
    vars,
    r_obs::NamedTuple,
    err_resp::MTResponse,
    mDist::MTModelDistribution,
    rDist::MTResponseDistribution;
    response_fields::Vector{Symbol}= [k for k ∈ fieldnames(typeof(rDist))],
    model_fields::Vector{Symbol}= [k for k ∈ fieldnames(typeof(mDist))],
    trans_utils::NamedTuple = (m = log_tf, h = lin_tf)
    ) where {model <: AbstractModel}

    m₀ = [];
    for k ∈ propertynames(mDist)
        if k in model_fields
            var_inf ~ getproperty(mDist, k)
            push!(m₀, broadcast(getproperty(trans_utils[k], :itf), var_inf));
        else
            push!(m₀, getfield(m_sample, k));   
        end
    end

    m_sample = typeof(m_sample)(m₀...);
    r_sample = forward(m_sample, vars);

    for k in response_fields
        r_obs[k] ~ getfield(rDist, k)(getfield(r_sample, k), getfield(err_resp, k))
    end
end