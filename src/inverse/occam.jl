"""
`occam_cache`: specifies the inverse algorithm while having a cache.
"""
mutable struct occam_cache{T}
    μgrid::Vector{T}
end
"""
`linsolve!`: Performs `inv(B)*y` using `LinearSolve.jl`
"""
function linsolve!(x, prob_init, B, y)
    prob_init.A = B
    prob_init.b = y
    x .= solve!(prob_init)
    nothing
end
"""
`Occam(;μgrid= [0.01, 1e6])`
"""
function Occam(; μgrid=[0.01, 1e6])
    return occam_cache{eltype(μgrid)}(μgrid)
end

"""
    function occam_step!(mₖ₊₁::model,
        respₖ₊₁::response,
        vars::Union{AbstractVector{Float32}, AbstractVector{Float64}},
        χ2::Union{Float64, Float32},
        μgrid::Vector{Float64},
        lin_utils::linear_utils,
        inv_utils::inverse_utils,
        trans_utils::transform_utils,
        linsolve_prob::LinearSolve.LinearCache;
        model_fields::Vector{Symbol}= [k for k ∈ fieldnames(typeof(mₖ₊₁))],
        response_fields::Vector{Symbol}= [k for k ∈ fieldnames(typeof(respₖ₊₁))],
        verbose= false
        ):

performs a single step of occam inversion, using golden line search.

### Variables:

  - `mₖ₊₁`: to store the next update, which will eventually be copied to mₖ
  - `respₖ₊₁`: to store the response for `mₖ₊₁`, for error calculation and anything
  - `vars`: to compute the forward model
  - `χ2::Union{Float64, Float32}`: threshold chi-squared error that needs to be met,
  - `μgrid::Vector{Float64}`: contains end points of the bounds for the lagrange multiplier,
  - `lin_utils::linear_utils`: contains the mₖ, Jₖ, Fₖ associate with the current iteration,
  - `inv_utils::inverse_utils`: contains D= ∂(n), W and dobs,
  - `trans_utils::transform_utils`: to  transform to and from the computational domain,
  - `linsolve_prob::LinearSolve.LinearCache`: for faster inverse operations,
  - `model_fields::Vector{Symbol`: which fields in model to consider changing,
  - `response_fields::Vector{Symbol}`: which fields in response to invert for, # to store the next update, which will eventually be copied to mₖ
  - `verbose`: whether to print the updates or not, default is true # model to regularize against
"""
function occam_step!(mₖ₊₁::model1, # to store the next update, which will eventually be copied to mₖ
        respₖ₊₁::response, # to store the response for mₖ₊₁, for error calculation and anything
        vars::Union{AbstractVector{Float32}, AbstractVector{Float64}}, # to compute the forward model
        χ2::Union{Float64, Float32}, # threshold chi-squared error that needs to be met
        μgrid::Vector{Float64}, # contains end points of the bounds for the lagrange multiplier
        lin_utils::linear_utils, # contains the mₖ, Jₖ, Fₖ associate with the current iteration
        inv_utils::inverse_utils, # contains D= ∂(n), W and dobs
        trans_utils::transform_utils, # to  transform to and from the computational domain
        linsolve_prob::LinearSolve.LinearCache; # for faster inverse operations
        model_fields::Vector{Symbol}=[k for k in fieldnames(typeof(mₖ₊₁))],
        response_fields::Vector{Symbol}=[k for k in fieldnames(typeof(respₖ₊₁))],
        verbose::Bool=true) where {
        model1 <: AbstractGeophyModel, model2 <: AbstractGeophyModel,
        response <: AbstractGeophyResponse}
    ϕ = (1 + sqrt(5)) / 2
    chi2min = (typeof(χ2))(1e6)
    μ = zero(eltype(μgrid))
    count = 0 # so that iterations do not run forever (will rarely happen, if it will)

    function f(x)
        linsolve!(mₖ₊₁.m,
            linsolve_prob,
            x .* inv_utils.D' * inv_utils.D .+ lin_utils.Jₖ' * inv_utils.W * lin_utils.Jₖ,
            lin_utils.Jₖ' *
            inv_utils.W *
            (inv_utils.dobs + lin_utils.Jₖ * lin_utils.mₖ - lin_utils.Fₖ))
        for k in model_fields # to model domain
            getfield(mₖ₊₁, k) .= 10.0 .^ trans_utils.tf.(getfield(mₖ₊₁, k))
        end
        forward!(respₖ₊₁, mₖ₊₁, vars)
        return χ²(reduce(vcat, [copy(getfield(respₖ₊₁, k)) for k in response_fields]),
            inv_utils.dobs; W=inv_utils.W)
    end

    x₁ = μgrid[1]
    x₃ = μgrid[end]
    x₂ = 10.0^((log10(x₃) + ϕ * log10(x₁)) / (1 + ϕ))
    x₄ = 10.0^((log10(x₁) + ϕ * log10(x₃)) / (1 + ϕ))

    # fx₁ = f(x₁)
    # fx₃ = f(x₃)
    fx₂ = f(x₂)
    fx₄ = f(x₄)

    tol = 1e-5
    count = 0
    while (x₃ - x₁) >= tol
        count += 1
        if count > 100
            verbose && (print("100 golden section iterations done. \t"))
            break
        end
        if fx₄ > fx₂
            x₃ = x₄
            x₄ = x₂
            fx₄ = fx₂
            x₂ = 10.0^((log10(x₃) + ϕ * log10(x₁)) / (1 + ϕ))
            fx₂ = f(x₂)

        else
            x₁ = x₂
            x₂ = x₄
            fx₂ = fx₄
            x₄ = 10.0^((log10(x₁) + ϕ * log10(x₃)) / (1 + ϕ))
            fx₄ = f(x₄)
        end
    end
    μ = sqrt(x₁ * x₃)

    # At the moment mₖ₊₁ contains the update for the last μ, we rewrite it with the best μ found.

    linsolve!(mₖ₊₁.m,
        linsolve_prob,
        μ .* inv_utils.D' * inv_utils.D .+ lin_utils.Jₖ' * inv_utils.W * lin_utils.Jₖ,
        lin_utils.Jₖ' *
        inv_utils.W *
        (inv_utils.dobs + lin_utils.Jₖ * lin_utils.mₖ - lin_utils.Fₖ))

    for k in model_fields # to model domain
        getfield(mₖ₊₁, k) .= 10.0 .^ trans_utils.tf.(getfield(mₖ₊₁, k)) # why do we have 10^ here
    end

    forward!(respₖ₊₁, mₖ₊₁, vars)

    verbose && (print("Works golden section search: μ= $μ, χ²= ",
        χ²(reduce(vcat, [copy(getfield(respₖ₊₁, k)) for k in response_fields]),
            inv_utils.dobs; W=inv_utils.W),
        "\n"))
    return μ
end

# we'd need a test sometime in future to check if the `r_obs` is indeed a response of `forward(m)`.
