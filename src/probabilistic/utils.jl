"""
    mutable struct struct mcmc_cache{T1 <: AbstractGeophyModelDistribution, T2 <: AbstractGeophyResponseDistribution}
        apriori::T1
        likelihood::T2
        n_samples::Int
        sampler
    end

placeholder to store

  - apriori in the form of any subtype of [`AbstractModelDistribution`](@ref)
  - likelihood in the form of any subtype of [`AbstractResponseDistribution`](@ref)
  - number of samples to obtain in `n_samples`
  - `Turing.jl` sampler to be used in `sampler`
"""
mutable struct mcmc_cache{
    T1 <: AbstractModelDistribution, T2 <: AbstractResponseDistribution}
    apriori::T1
    likelihood::T2
    n_samples::Int
    sampler
end

"""
    @model function mcmc_turing(
        m_sample::model,
        vars,
        r_obs::NamedTuple,
        err_resp::MTResponse,
        mDist::mdist,
        rDist::rdist;
        response_fields::Vector{Symbol}= [k for k ∈ fieldnames(typeof(rDist))],
        model_fields::Vector{Symbol}= [k for k ∈ fieldnames(typeof(mDist))],
        trans_utils::NamedTuple = (m = log_tf, h = lin_tf)
        ) where {model <: AbstractModel, mdist <: AbstractModelDistribution, rdist <: AbstractResponseDistribution}

makes a `Turing.jl` model to perform MCMC sampling

### Variables:

  - `vars`: variables that need to be passed into the `forward` function along with `model` to generate a `response`
  - `r_obs`: named tuple containing the observed data, with the same keys as the fields in the corresponding `response`
  - `err_resp`: `response` variable that contains the errors
  - `mDist`: any subtype of [`AbstractModelDistribution`](@ref) contains the apriori information
  - `rDist`: any subtype of [`AbstractResponseDistribution`](@ref) contains the likelihood information

### Keyword/optional arguments

  - `response_fields`:  which fields in `response` to invert for
  - `model_fields`: fields in `model` to draw inference on
  - `trans_utils`: to transform the model field variables to and from computational (inference) domain
"""
@model function mcmc_turing(m_sample::model,
        const_data,
        vars,
        r_obs::NamedTuple,
        err_resp::response,
        mDist::mdist,
        rDist::rdist;
        response_fields::Vector{Symbol}=[k for k in fieldnames(typeof(rDist))],
        model_fields::Vector{Symbol}=[k for k in fieldnames(typeof(mDist))],
        model_trans_utils::NamedTuple=(m=lin_tf, h=lin_tf),
        response_trans_utils::NamedTuple=(ρₐ=lin_tf, ϕ=lin_tf)) where {
        model <: AbstractModel, response <: AbstractResponse,
        mdist <: AbstractModelDistribution, rdist <: AbstractResponseDistribution}
    m0 = (; zip([propertynames(mDist)...], const_data)...)

    for k in propertynames(mDist)
        if k in model_fields
            m0[k] ~ getproperty(mDist, k)
        end
    end

    m_sample = typeof(m_sample)([broadcast(getproperty(model_trans_utils[k], :tf), m0[k])
                                 for k in propertynames(mDist)]...)

    r_sample = forward(m_sample, vars; trans_utils=response_trans_utils)

    for k in response_fields
        r_obs[k] ~ getfield(rDist, k)(getfield(r_sample, k), getfield(err_resp, k) .^ 2)
    end
end
