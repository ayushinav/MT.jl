"""
`nl_cache`: specifies the inverse algorithm while having a cache.
"""
mutable struct nl_cache{T1, T2}
    alg::T1
    μ::T2
    kwargs::NamedTuple
end

mutable struct NonlinearAlg{T1, T2}
    alg::T1
    μ::T2
    kwargs::NamedTuple
end

# ===== NonlinearSolve.jl =========

"""
    function inverse!(mₖ::model,
            robs::response,
            vars::Vector{Float64},
            alg_cache::nl_cache;
            W= nothing, # Weight matrix
            max_iters= 20, χ2=1.,
            response_fields::Vector{Symbol}= [k for k ∈ fieldnames(typeof(robs))],
            model_fields::Vector{Symbol}= [k for k ∈ fieldnames(typeof(m₀))], # this will not be used but for the sake of generality for all inverse algs
            trans_utils::transform_utils= log_sigmoid_tf,
            verbose= true
        ):

updates `mₖ` using occam iteration to fit `robs` within a misfit of `χ2`, by default set to 1.0.

### Variables:

  - `mₖ::model`: Initial model guess, will be updated during the inverse process
  - `robs::response`: response to invert for
  - `vars::Vector{Float64}`: variables required for forward modeling, eg., `ω` for MT
  - `alg_cache::occam_cache`: deterimines the algorithm to be performed for inversion
  - `W= nothing`: Weight matrix, will be `I` if nothing is provided
  - `max_iters= 20`: maximum number of iterations
  - `χ2=1.`: threshold misfit
  - `response_fields::Vector{Symbol}= [k for k ∈ fieldnames(typeof(robs))]`: choose data of response to perform inversion on, eg., ρₐ for MT, by default chooses all the data (ρₐ and ϕ)
  - `model_fields::Vector{Symbol}= [k for k ∈ fieldnames(typeof(m₀))]`: will generally be fixed, see docs for details
  - `trans_utils::transform_utils= sigmoid_tf`: for bounds transformation,
  - `verbose`: whether to print updates after each iteration, defaults to true

### Returns:

return message in the form of `return_code` and updates `mₖ` in-place.

### Example:

`inverse!(m_occam, r_obs, Occam([1e-2, 1e6]))`
"""
function inverse!(mₖ::model1,
        robs::response,
        vars::Vector{Float64},
        alg_cache::NonlinearAlg;
        W=nothing,
        L=nothing,
        max_iters=30,
        χ2=1.0,
        response_fields::Vector{Symbol}=[k for k in fieldnames(typeof(robs))],
        # model_fields::Vector{Symbol}=[k for k in fieldnames(typeof(mₖ))], # this will not be used but for the sake of generality for all inverse algs
        model_trans_utils::transform_utils=pow_sigmoid_tf,
        response_trans_utils::transform_utils=log_tf,
        verbose::Bool=true,
        mᵣ=nothing) where {
        model1 <: AbstractGeophyModel, response <: AbstractGeophyResponse}
    prec = eltype(mₖ.m)
    model_fields = [:m]

    (W === nothing) && (W = prec.(I(n_resp)))
    (W === nothing) && (W = prec.(∂(length(mₖ.m))))

    model_type = typeof(mₖ).name.wrapper

    p = (
        model_type = model_type, 
        h = mₖ.h, 
        model_trans_utils = model_trans_utils, 
        vars = vars, 
        response_fields = response_fields, W = W, μ = alg_cache.μ,
        r_obs = robs, L = L
    )

    # construct_cost_function(mₖ.m, p)

    prob = SciMLBase.NonlinearLeastSquaresProblem(construct_cost_function, model_trans_utils.itf.(mₖ.m), p)
    nlcache = init(prob, alg_cache.alg(); alg_cache.kwargs...)
    iters = 1
    chi2 = 1e6

    @show nlcache.u

    while iters <= max_iters

        # check misfit condition
        model = model_type(model_trans_utils.tf.(nlcache.u), mₖ.h)
        resp_ = forward(model, vars)
        chi2 = χ²(reduce(vcat, [getfield(resp_, k) for k in response_fields]),
        reduce(vcat, [getfield(robs, k) for k in response_fields]); W=W)
        @show iters, chi2
        if chi2 <= χ2
            break;
        end

        step!(nlcache)
        # @show nlcache.u

        iters += 1 
    end

    mₖ.m .= model_trans_utils.tf.(nlcache.u)

    return return_code(chi2 <= χ2, (μ=alg_cache.μ,), mₖ, χ2, chi2)
end

# focussing on geophysical models for now
# Not performant at the moment
function construct_cost_function(m, p) 
    @unpack model_type, h, model_trans_utils, vars, response_fields, W, μ, r_obs, L = p
    # model = model_type(model_trans_utils.tf.(m), h)
    model = model_type(broadcast(model_trans_utils.tf, m), h)
    resp_ = forward(model, vars)

    L1 = χ²(reduce(vcat, [getfield(resp_, k) for k in response_fields]),
    reduce(vcat, [getfield(r_obs, k) for k in response_fields]); W=W)
    L2 = μ * sqrt(norm(L * m))
    @show L1, L2

    return [L1 + L2]
end
