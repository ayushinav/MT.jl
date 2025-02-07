
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
        trans_utils::NamedTuple=(m=sigmoid_tf, h=lin_tf)) where {
        resp1 <: AbstractGeophyResponse, resp2 <: AbstractGeophyResponse}

    # first an occam iteration to get the first model
    W = Diagonal(vcat([inv.(getfield(err_resp, k)) .^ 2
                       for k in fieldnames(typeof(err_resp))]...))

    # retcode = inverse!(alg_cache.m₀, r_obs, vars, occam_cache(alg_cache.μgrid);
    #     W=W, χ2=alg_cache.χ2, max_iters=alg_cache.max_iters_occam,
    #     response_fields=alg_cache.response_fields, verbose=alg_cache.verbose)

    # @show retcode.misfit_achieved

    pert_model = copy(alg_cache.m₀)

    pert_resp = copy(r_obs)

    n = length(alg_cache.m₀.m)
    L = ∂(n)

    μ = 1.0
    resp_ = copy(r_obs)
    ϕ = (1 + sqrt(5)) / 2

    m_chains = zeros(n, alg_cache.n_samples)
    μ_chains = zeros(1, alg_cache.n_samples)

    for i in 1:(alg_cache.n_samples)
        (i % 5 == 0) && (@show i)

        ## Step 1

        # perturbed response
        for k in fieldnames(typeof(r_obs))
            getfield(pert_resp, k) .= getfield(r_obs, k) .+
                                      randn(size(getfield(err_resp, k))) .*
                                      getfield(err_resp, k)
        end

        # perturbed model : the one we regularize against
        # we first draw a perturbation in the computational domain

         #reg term
        mul!(pert_model.m, L, randn(eltype(pert_model.m), size(pert_model.m))) # L * ξ
        rmul!(pert_model.m, sqrt(μ)) # sqrt(μ) * L * ξ 

        # pulling everything from computational domain

        # broadcast!((x) -> (10.0^trans_utils[:m].tf(x)), pert_model.m, pert_model.m) # moving to model domain, we put everything to computational domain in inv.jl
        # 10.0 .^ trans_utils.tf.(getfield(mₖ₊₁, k))
        fill!(alg_cache.m₀.m, 100.0)

        # @show alg_cache.m₀.m
        # @show pert_model.m
        ret_code = inverse!(
            alg_cache.m₀, r_obs, vars, occam_cache([μ, μ]); W=W, χ2=alg_cache.χ2,
            max_iters=alg_cache.max_iters_occam, response_fields=alg_cache.response_fields,
            verbose=alg_cache.verbose, mᵣ=pert_model, trans_utils=trans_utils[:m])

        m_chains[:, i] .= alg_cache.m₀.m

        # @show μ
        # @show alg_cache.m₀.m

        ## Step 2

        # perturbed response
        for k in fieldnames(typeof(r_obs))
            getfield(pert_resp, k) .= getfield(r_obs, k) .+
                                      randn(size(getfield(err_resp, k))) .*
                                      getfield(err_resp, k)
        end

        broadcast!((x) -> (trans_utils[:m].itf(log10.(x))), alg_cache.m₀.m, alg_cache.m₀.m) # to computational domain
        rmul!(alg_cache.m₀.m, sqrt(μ)) # current μ

        function g(x) #, alg_cache, trans_utils, resp_, vars)
            rmul!(alg_cache.m₀.m, sqrt(inv(x)))

            # broadcast!(trans_utils[:m].itf, alg_cache.m₀.m, alg_cache.m₀.m) # to model domain
            # broadcast!((x) -> (10.0^trans_utils[:m].tf(x)), alg_cache.m₀.m, alg_cache.m₀.m) #  to model domain
            alg_cache.m₀.m .= 10 .^ trans_utils[:m].tf.(alg_cache.m₀.m)

            forward!(resp_, alg_cache.m₀, vars)
            # @show "chi2_err"
            chi2_err = χ²(
                reduce(vcat, [getfield(resp_, k) for k in alg_cache.response_fields]),
                reduce(vcat, [getfield(r_obs, k) for k in alg_cache.response_fields]);
                W=W)

            # broadcast!(trans_utils[:m].tf, alg_cache.m₀.m, alg_cache.m₀.m) # to computational domain
            alg_cache.m₀.m .= trans_utils[:m].itf.(log10.(alg_cache.m₀.m))
            # broadcast!(
            #     (x) -> (trans_utils[:m].itf(log10.(x))), alg_cache.m₀.m, alg_cache.m₀.m) # to computational domain
            rmul!(alg_cache.m₀.m, sqrt(x))

            return chi2_err
        end

        x₁ = alg_cache.μgrid[1]
        x₃ = alg_cache.μgrid[end]
        x₂ = 10.0^((log10(x₃) + ϕ * log10(x₁)) / (1 + ϕ))
        x₄ = 10.0^((log10(x₁) + ϕ * log10(x₃)) / (1 + ϕ))

        # fx₁ = f(x₁)
        # fx₃ = f(x₃)
        fx₂ = g(x₂)
        fx₄ = g(x₄)

        tol = 1e-5
        count = 0
        while (x₃ - x₁) >= tol
            count += 1
            if count > 100
                print("100 golden section iterations done. \t")
                break
            end
            if fx₄ > fx₂
                x₃ = x₄
                x₄ = x₂
                fx₄ = fx₂
                x₂ = 10.0^((log10(x₃) + ϕ * log10(x₁)) / (1 + ϕ))
                fx₂ = g(x₂)

            else
                x₁ = x₂
                x₂ = x₄
                fx₂ = fx₄
                x₄ = 10.0^((log10(x₁) + ϕ * log10(x₃)) / (1 + ϕ))
                fx₄ = g(x₄)
            end
        end
        μ = sqrt(x₁ * x₃)

        # copyto!(alg_cache.m₀.m, m_chains[:, i])
        μ_chains[1, i] = μ

        # @show μ
        # @show m_chains[:, i]
    end

    return Turing.Chains(vcat(m_chains, μ_chains)', [Symbol("m[$i]") for i in 1:(n + 1)])
end

# mutable struct RTO_MTModel <: AbstractGeophyModel
#     ξ
#     h
#     μ
# end
