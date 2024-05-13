@testset "inverse tests" begin
    h= [1000., 1000.]; # m
    ρ= [100., 10., 1000.]; # Ωm
    m= MTModel(ρ, h);
    T= 10 .^(range(-2, 2, length= 57));
    ω= 2π./T;
    nω= length(T);

    r_obs = forward(m, ω);
    m_occam = MTModel(fill(100., 50), collect(range(0, 5e3, length = 49)))


    err_ρ= 0.1 .* r_obs.ρₐ
    err_ϕ= asin(0.01)* 180/π .* ones(length(ω))
    err_resp= diagm([err_ρ..., err_ϕ...])

    ret_code = inverse!(m_occam, r_obs, ω, Occam(μgrid = [1e-2, 1e6]), W = err_resp, χ2 = 1., max_iters = 50, verbose = false);

    @test ret_code.if_pass === true
    @test ret_code.if_pass == true
end