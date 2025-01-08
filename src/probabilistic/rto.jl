
struct rto_cache{model <: AbstractGeophyModel, T1 <: AbstractVector{<:Any},
    T2 <: Union{Float64, Float32}}
    m₀::model
    μgrid::T1
    max_iters_occam::Int
    n_samples::Int
    χ2::T2
    response_fields::Vector{Symbol}
    verbose::Bool
end

function stochastic_inverse(r_obs::resp1,
        err_resp::resp2,
        vars,
        alg_cache::rto_cache;
        trans_utils::NamedTuple=(m=log_tf, h=lin_tf)) where {
        resp1 <: AbstractGeophyResponse, resp2 <: AbstractGeophyResponse}

    # first an occam iteration to get the first model
    W = Diagonal(vcat([inv.(getfield(err_resp, k)) .^ 2
                       for k in fieldnames(typeof(err_resp))]...))

    retcode = inverse!(alg_cache.m₀, r_obs, vars, occam_cache(alg_cache.μgrid);
        W=W, χ2=alg_cache.χ2, max_iters=alg_cache.max_iters_occam,
        response_fields=alg_cache.response_fields, verbose=alg_cache.verbose)

    # alg_cache.m₀ is in model domain, not computational domain

    @show retcode.misfit_achieved
    pert_model = copy(alg_cache.m₀)
    pert_m = zero(alg_cache.m₀.m)

    pert_resp = copy(r_obs)

    n = length(alg_cache.m₀.m)
    L = ∂(n)

    # check this L'L
    lin_prob = LinearProblem(L'L, pert_m)
    linsolve_prob = init(lin_prob; assumptions=LinearSolve.OperatorAssumptions(true)) #; condition=LinearSolve.OperatorCondition.WellConditioned))

    μ = 1.0

    m_chains = zeros(n, alg_cache.n_samples)

    # @show alg_cache.m₀.m

    @showprogress for i in 1:(alg_cache.n_samples)

        ## RTO step

        # perturbed response
        for k in fieldnames(typeof(r_obs))
            getfield(pert_resp, k) .= rand(MultivariateNormal(
                getfield(r_obs, k), Diagonal(getfield(err_resp, k))))
        end

        # perturbed model : the one we regularize against
        # we first draw a perturbation in the computational domain

        pert_m .= rand(MultivariateNormal(
            zero(pert_m), Diagonal(ones(n) .|> eltype(pert_m)))) # sample ξ

        rmul!(pert_m, inv(sqrt(μ)))
        linsolve!(pert_model.m, linsolve_prob, L, pert_m)

        # next steps are strictly making things geophysical here
        # pulling everything from computational domain

        broadcast!(trans_utils[:m].tf, pert_model.m, pert_model.m) # we put everything to computational domain in inv.jl

        ret_code = inverse!(alg_cache.m₀, r_obs, vars, occam_cache([μ, μ]); W=W,
            χ2=alg_cache.χ2, max_iters=alg_cache.max_iters_occam,
            response_fields=alg_cache.response_fields,
            verbose=alg_cache.verbose, mᵣ=pert_model)

        ## TKO step

        # perturbed response
        for k in fieldnames(typeof(r_obs))
            getfield(pert_resp, k) .= rand(MultivariateNormal(
                getfield(r_obs, k), Diagonal(getfield(err_resp, k))))
        end

        function f(x) end

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

        # the next occam function needs m to be a function of μ and DOES NOT INCLUDE ANY REGULARIZATION
        # perturbed response
        for k in fieldnames(typeof(r_obs))
            getfield(pert_resp, k) .= rand(MultivariateNormal(
                getfield(r_obs, k), Diagonal(getfield(err_resp, k))))
        end

        ret_code = inverse!(copy(alg_cache.m₀), r_obs, vars, # 23 use copy(m₀) because we don't want to perturb m₀
            occam_cache(alg_cache.μgrid);
            W=W, χ2=alg_cache.χ2, max_iters=1,
            response_fields=alg_cache.response_fields, verbose=alg_cache.verbose)

        μ = ret_code.parameters[:μ]

        m_chains[:, i] .= alg_cache.m₀.m
    end

    return Turing.Chains(
        broadcast(trans_utils[:m].itf, m_chains)', [Symbol("ρ[$i]") for i in 1:n])
end

# mutable struct RTO_MTModel <: AbstractGeophyModel
#     ξ
#     h
#     μ
# end

#= RTO issue
We know about observed data d_obs, the errors associated in C\_

We first use Occam to get m⁰

=#
