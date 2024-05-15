
struct rto_cache{model<:AbstractGeophyModel, T1<:AbstractVector{<:Any}, T2<:Union{Float64, Float32}}
    m₀::model
    μgrid::T1
    max_iters_occam::Int
    n_samples::Int
    χ2::T2
    response_fields::Vector{Symbol}
    verbose::Bool

end

function stochastic_inverse(
    r_obs::resp1,
    err_resp::resp2,
    vars,
    alg_cache::rto_cache;
    trans_utils::NamedTuple = (m = logg_tf, h = lin_tf)
    ) where {resp1 <: AbstractGeophyResponse, resp2 <: AbstractGeophyResponse}


    # first an occam iteration to get the first model
    W = Diagonal(
        vcat([inv.(getfield(err_resp, k)).^2 for k in fieldnames(typeof(err_resp))]...)
    );

    retcode = inverse!(alg_cache.m₀, r_obs, vars, 
        occam_cache(alg_cache.μgrid), W = W, χ2 = alg_cache.χ2,
        max_iters = alg_cache.max_iters_occam,
        response_fields = alg_cache.response_fields, 
        verbose = alg_cache.verbose
    )

    @show retcode.if_pass
    pert_model = copy(alg_cache.m₀);
    pert_m = zero(alg_cache.m₀.m);

    pert_resp = copy(r_obs);

    n= length(alg_cache.m₀.m);
    L = MT.∂(n);
    lin_prob = LinearProblem(L, pert_m);
    linsolve_prob = init(lin_prob, assumptions=LinearSolve.OperatorAssumptions(true, condition= LinearSolve.OperatorCondition.WellConditioned));

    μ = 1.0;

    m_chains = zeros(n, alg_cache.n_samples);

    # @show alg_cache.m₀.m

    @showprogress for i in 1:alg_cache.n_samples
        
       ## RTO step

       # perturbed response
        for k in fieldnames(typeof(r_obs))
            getfield(pert_resp, k) .= rand(
                MultivariateNormal(getfield(r_obs, k), Diagonal(getfield(err_resp, k)))
            )
        end

        # perturbed model : the one we regularize against
        # we first draw a pertubation in the computational domain

        pert_m .= rand(MultivariateNormal(zero(pert_m), 
            Diagonal(ones(n) .|> eltype(pert_m))
        ))

        rmul!(pert_m, inv(sqrt(μ)));
        linsolve!(pert_model.m, linsolve_prob, L, pert_m);

        # now this is strictly making things geophysical here.
        # pulling everything from computational domain

        broadcast!(trans_utils[:m].tf, pert_model.m, pert_model.m);

        ret_code = inverse!(
            alg_cache.m₀, r_obs, vars, 
            occam_cache([μ, μ]), W = W, χ2 = alg_cache.χ2,
            max_iters = alg_cache.max_iters_occam,
            response_fields = alg_cache.response_fields, 
            verbose = alg_cache.verbose, 
            mᵣ = pert_model
        );

        ## TKO step

        # perturbed response
        for k in fieldnames(typeof(r_obs))
            getfield(pert_resp, k) .= rand(
                MultivariateNormal(getfield(r_obs, k), Diagonal(getfield(err_resp, k)))
            )
        end

        ret_code = inverse!(
            copy(alg_cache.m₀), r_obs, vars, # 23 use copy(m₀) because we don't want to perturb m₀
            occam_cache(alg_cache.μgrid), W = W, χ2 = alg_cache.χ2,
            max_iters = 1,
            response_fields = alg_cache.response_fields, 
            verbose = alg_cache.verbose 
        );

        μ = ret_code.parameters[:μ];

        m_chains[:, i].= alg_cache.m₀.m;

    end

    return Turing.Chains(broadcast(trans_utils[:m].itf, m_chains)', [Symbol("ρ[$i]") for i in 1:n]);

end