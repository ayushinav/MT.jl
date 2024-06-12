using Distributions
using Turing

@testset "mcmc_tests" begin

    m_test = MTModel([100., 10., 1000.], [1e3, 1e3]);
    f = 10 .^ range(-4, stop = 1, length = 25);
    ω = vec(2π .* f);

    r_obs = forward(m_test, ω);

    err_phi = asin(0.01) * 180/π .* ones(length(ω));
    err_appres = 0.02 * r_obs.ρₐ;
    err_resp = MTResponse(err_appres, err_phi);

    r_obs.ρₐ .= r_obs.ρₐ .+ err_appres;
    r_obs.ϕ .= r_obs.ϕ .+ err_phi;

    respD = MTResponseDistribution(normal_dist, normal_dist);

    z = 10 .^collect(range(1, stop = 4, length = 100));
    h = diff(z);


    modelD = MTModelDistribution(
        Product(
        [Uniform(-1., 5.) for i in eachindex(z)]
        ),
        vec(h)
    );


    n_samples = 50;
    mcache = mcmc_cache(modelD, respD, 50, NUTS());

    log_tf2= transform_utils([], log10, (x) -> 10^x, (x) -> inv(x * log(10)));

    mcmc_chain = stochastic_inverse(r_obs, err_resp, ω, mcache, trans_utils = (m = log_tf,))

    preds= []
    for k in mcmc_chain.name_map.parameters
        push!(preds, mcmc_chain[k].data[:])
    end
    pred= hcat(preds...);
    pred.= 10 .^ pred;

    rho_err = 0.;
    ph_err = 0.;
    for idx in 1:n_samples
        rand_model = pred[idx, :];

        m1_model = rand_model[:];
        h_model = copy(h); 

        m_model = MTModel(m1_model, h_model);
        resp_model = forward(m_model, ω);

        rho_err += sum(abs.(resp_model.ρₐ .- r_obs.ρₐ)) 
        ph_err += sum(abs.(resp_model.ϕ .- r_obs.ϕ)) 
    end

    @test sqrt(rho_err)/n_samples .<= sum(err_resp.ρₐ)
    @test sqrt(ph_err)/n_samples .<= sum(err_resp.ϕ)
end