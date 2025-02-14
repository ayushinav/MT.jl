"""
    function inverse!(mₖ::model,
            robs::response,
            vars::Vector{Float64},
            alg_cache::occam_cache;
            W= nothing,
            L= nothing,
            max_iters= 30, χ2=1.,
            response_fields::Vector{Symbol}= [k for k ∈ fieldnames(typeof(robs))],
            model_fields::Vector{Symbol}= [k for k ∈ fieldnames(typeof(m₀))], # this will not be used but for the sake of generality for all inverse algs
            model_trans_utils::transform_utils= sigmoid_tf,
            verbose= true
        ):

updates `mₖ` using occam iteration to fit `robs` within a misfit of `χ2`, by default set to 1.0

### Variables:

  - `mₖ`: Initial model guess, will be updated during the inverse process
  - `robs`: response to invert for
  - `vars`: variables required for forward modeling, eg., `ω` for MT
  - `alg_cache`: deterimines the algorithm to be performed for inversion
  - `W`: Weight matrix, defaults to identity matrix `I`
  - `L`: Regularization matrix, defaults to derivative matrix, given by `∂`(@ref)
  - `max_iters= 30`: maximum number of iterations
  - `χ2=1.`: target misfit
  - `response_fields: choose data of response to perform inversion on, eg., ρₐ for MT, by default chooses all the data (ρₐ and ϕ)
  - `model_trans_utils:transform_utils= pow_sigmoid_tf`: conversion to and from computational domain,
  - `response_trans_utils`: for scaling the response parameters,
  - `mᵣ`: model in physical domain to be regularized against
  - `reg_term`: (For internals) When model in physical domain does not exist, `reg_term` helps, eg. case of RTO-TKO
  - `verbose`: whether to print updates after each iteration, defaults to true

### Returns:

return message in the form of `return_code` and updates `mₖ` in-place.

### Example:

`inverse!(m_occam, r_obs, Occam([1e-2, 1e6]))`
"""
function inverse!(mₖ::model1,
        robs::response,
        vars::Vector{Float64},
        alg_cache::occam_cache;
        W=nothing,
        L=nothing,
        max_iters=30,
        χ2=1.0,
        response_fields::Vector{Symbol}=[k for k in fieldnames(typeof(robs))],
        model_trans_utils::transform_utils=sigmoid_tf,
        response_trans_utils::NamedTuple=(; ρₐ=lin_tf, ϕ=lin_tf),
        mᵣ=nothing,
        reg_term=nothing,
        verbose::Bool=true) where {
        model1 <: AbstractGeophyModel, response <: AbstractGeophyResponse}
    prec = eltype(mₖ.m)
    model_fields = [:m]

    n_model = length(mₖ.m)
    n_vars = length(vars)
    n_resp = length(response_fields) * n_vars

    (W === nothing) && (W = prec.(I(n_resp)))
    (L === nothing) && (L = prec.(∂(n_model)))

    lin_utils = linear_utils(
        view(mₖ.m, :), zeros(prec, n_resp), zeros(prec, n_resp, n_model))

    respₖ = zero_abstract(robs)
    jc = jacobian_cache(response_fields, robs, mₖ, model_fields)

    for (i, k) in enumerate(response_fields)
        setfield!(jc.j, k, view(lin_utils.Jₖ, ((i - 1) * n_vars + 1):(i * n_vars), :))
        setfield!(respₖ, k, view(lin_utils.Fₖ, ((i - 1) * n_vars + 1):(i * n_vars)))
    end

    inv_utils = inverse_utils(
        L, W, reduce(vcat, [copy(getfield(robs, k)) for k in response_fields]))

    mₖ₊₁ = copy(mₖ)
    respₖ₊₁ = copy(respₖ)

    lin_prob = LinearProblem(inv_utils.D'inv_utils.D,
        lin_utils.Jₖ' * (inv_utils.dobs + lin_utils.Jₖ * lin_utils.mₖ - lin_utils.Fₖ))
    linsolve_prob = init(lin_prob;
        assumptions=LinearSolve.OperatorAssumptions(
            true; condition=LinearSolve.OperatorCondition.WellConditioned))

    forward!(respₖ, mₖ, vars; response_trans_utils=response_trans_utils) # for the first iteration
    itr = 1
    chi2 = prec(1e6)

    if mᵣ !== nothing
        for k in model_fields # to computational domain
            getfield(mᵣ, k) .= model_trans_utils.itf.((getfield(mᵣ, k)))
        end
    end

    if reg_term === nothing
        reg_term = zero(mₖ.m)
    end

    μ_last = 0.0
    while itr <= max_iters
        verbose && (print("$itr: "))
        jacobian!(mₖ, vars, jc; model_fields=model_fields, response_fields=response_fields)

        for k in model_fields # to computational domain
            getfield(mₖ, k) .= model_trans_utils.itf.(getfield(mₖ, k))
        end

        μ_last = occam_step!(mₖ₊₁, # to store the next update, which will eventually be copied to mₖ
            respₖ₊₁, # to store the response for mₖ₊₁, for error calculation and anything
            vars, # to compute the forward model
            χ2, # threshold chi-squared error that needs to be met
            alg_cache.μgrid, # for gridsearch of μ for Occam
            lin_utils, # contains the mₖ, Jₖ, Fₖ associate with the current iteration
            inv_utils, # contains D= ∂(n), W and dobs
            model_trans_utils, # to  transform to and from the computational domain
            response_trans_utils, linsolve_prob; # for faster inverse operations
            model_fields=model_fields, response_fields=response_fields,
            mᵣ=mᵣ, verbose=verbose, reg_term=reg_term)

        for k in model_fields # copying things to mₖ
            getfield(mₖ, k) .= getfield(mₖ₊₁, k)
        end
        forward!(respₖ, mₖ, vars; response_trans_utils=response_trans_utils)
        chi2 = χ²(reduce(vcat, [copy(getfield(respₖ, k)) for k in response_fields]),
            inv_utils.dobs; W=inv_utils.W)
        if chi2 < χ2
            break
        end
        itr += 1
    end

    if mᵣ !== nothing
        for k in model_fields # back to model domain
            getfield(mᵣ, k) .= model_trans_utils.tf.((getfield(mᵣ, k)))
        end
    end

    return return_code(chi2 <= χ2, (μ=μ_last,), mₖ, χ2, chi2)
end
