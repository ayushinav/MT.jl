@testitem "conductivity tests" tags = [:rp] begin
    using MT
    methods_list = [
        SEO3,
        UHO2014,
        Jones2012,
        Yoshino2009,
        Wang2006,
        Poe2010,

        # melt
        Ni2011,
        Sifre2014,
        Gaillard2008,
    ]

    T = 1273.0f0:30:1573.0f0
    Ch2o_ol = 2.0f4:2.0f3:4.0f4
    Ch2o_m = 2.0f4:2.0f3:4.0f4
    Cco2_m = 1.0f4:2.0f3:3.0f4

    inps = (;
        SEO3 = [T],
        UHO2014 = [T, Ch2o_ol],
        Jones2012 = [T, Ch2o_ol],
        Yoshino2009 = [T, Ch2o_ol],
        Wang2006 = [T, Ch2o_ol],
        Poe2010 = [T, Ch2o_ol],

        # melt
        Ni2011 = [T, Ch2o_m],
        Sifre2014 = [T, Ch2o_m, Cco2_m],
        Gaillard2008 = [T],
    )

    outs = (;
        SEO3=Float32.([-4.0572, -3.8956, -3.7387, -3.5857, -3.4355, -3.2874, -3.1402,
            -2.9933, -2.8457, -2.6968, -2.5463]),
        UHO2014=Float32.([1.2728, 1.4152, 1.5459, 1.6664, 1.7780, 1.8818,
            1.9786, 2.0692, 2.1542, 2.2341, 2.3095]),
        Jones2012=Float32.([0.1550, 0.2774, 0.3920, 0.4996, 0.6011, 0.6971,
            0.7880, 0.8744, 0.9567, 1.0351, 1.1099]),
        Yoshino2009=Float32.([-0.6429, -0.5108, -0.3878, -0.2729, -0.1650,
            -0.0634, 0.0325, 0.1232, 0.2093, 0.2912, 0.3691]),
        Wang2006=Float32.([-0.3832, -0.2753, -0.1734, -0.0768, 0.0150,
            0.1023, 0.1857, 0.2653, 0.3415, 0.4144, 0.4844]),
        Poe2010=Float32.([3.4473, 3.6441, 3.8203, 3.9789, 4.1224, 4.2527,
            4.3714, 4.4800, 4.5796, 4.6711, 4.7554]),
        Ni2011=Float32.([-2.3579, -1.3975, -0.7500, -0.2847, 0.0652, 0.3375,
            0.5552, 0.7329, 0.8807, 1.0053, 1.1117]),
        Sifre2014=Float32.([-0.0128, 0.1570, 0.3119, 0.4537, 0.5837, 0.7031,
            0.8131, 0.9145, 1.0084, 1.0954, 1.1763]),
        Gaillard2008=Float32.([2.2276, 2.2577, 2.2865, 2.3140, 2.3403, 2.3655,
            2.3897, 2.4129, 2.4352, 2.4566, 2.4772])
    )

    @testset "$(methods_list[i])" for i in eachindex(methods_list)
        m = keys(inps)[i]
        model = methods_list[i](inps[m]...)
        out_ = forward(model, [])
        @inferred forward(model, [])
        @report_call forward(model, [])
        @test all(isapprox.(out_.σ, outs[m], atol=1e-2))
    end

    # mixing models
    # mixing_list = [HS1962_plus(), HS1962_minus(), MAL(0.2)]
    # mixing_outs = log10.([3.869f-4, 1.1502f-4, 2.8f-3])

    # @testset "$(mixing_list[i])" for i in eachindex(mixing_list)
    #     model = construct_mixing_models([1000.0 + 273.0, 2e4],
    #         [:T, :Ch2o_m],
    #         [0.1],
    #         [SEO3, Ni2011],
    #         [mixing_list[i]])

    #     out_ = forward(model, [])
    #     @test round(first(out_.σ); digits=2) ≈ round(mixing_outs[i]; digits=2)
    # end
end

@testitem "elasticity tests" tags = [:rp] begin
    T = collect(1273.0f0:30:1573.0f0)
    ρ = collect(3300.0f0:100.0f0:4300.0f0)
    ϕ = collect(1.0f-2:1.0f-3:2.0f-2)
    P = 2 .+ zero(T)

    inps = (
        anharmonic=[T, P, ρ],
        anharmonic_poro=[T, P, ρ, ϕ],
        SLB2005=[T, P]
    )

    outs = (
        anharmonic=MT.RockphyElastic(
            [7.1367f10, 7.0959f10, 7.0551f10, 7.0143f10, 6.9735f10, 6.9327f10,
                6.8919f10, 6.8511f10, 6.8103f10, 6.7695f10, 6.7287f10],
            [1.1895f11, 1.1827f11, 1.1759f11, 1.1691f11, 1.1623f11, 1.1555f11,
                1.1487f11, 1.1419f11, 1.1351f11, 1.1283f11, 1.1215f11],
            [8.0548f03, 7.9127f03, 7.7764f03, 7.6454f03, 7.5194f03, 7.3981f03,
                7.2811f03, 7.1682f03, 7.0591f03, 6.9537f03, 6.8516f03],
            [4.6504f03, 4.5684f03, 4.4897f03, 4.4141f03, 4.3413f03, 4.2713f03,
                4.2038f03, 4.1386f03, 4.0756f03, 4.0147f03, 3.9558f03]
        ), anharmonic_poro=RockphyElastic(
            [6.9053f10, 6.8464f10, 6.7877f10, 6.7293f10, 6.6711f10,
                6.6132f10, 6.5556f10, 6.4982f10, 6.4410f10, 6.3840f10, 6.3272f10],
            [1.1682f11, 1.1596f11, 1.1510f11, 1.1425f11, 1.1341f11,
                1.1257f11, 1.1173f11, 1.1090f11, 1.1007f11, 1.0925f11, 1.0843f11],
            [7.9531f03, 7.8040f03, 7.6610f03, 7.5236f03, 7.3914f03,
                7.2641f03, 7.1414f03, 7.0230f03, 6.9086f03, 6.7980f03, 6.6910f03],
            [4.5744f03, 4.4874f03, 4.4038f03, 4.3235f03, 4.2462f03,
                4.1717f03, 4.0999f03, 4.0306f03, 3.9635f03, 3.8987f03, 3.8359f03]
        ), SLB2005=RockphyElastic(
            zeros(Float32, 11),
            zeros(Float32, 11),
            zeros(Float32, 11),
            [4.4782f03, 4.4669f03, 4.4555f03, 4.4442f03, 4.4328f03, 4.4215f03,
                4.4102f03, 4.3988f03, 4.3875f03, 4.3761f03, 4.3648f03]
        )
    )

    methods_list = [anharmonic, anharmonic_poro, SLB2005]

    for i in eachindex(methods_list)
        m = keys(inps)[i]
        model = methods_list[i](inps[m]...)
        out_ = forward(model, [])
        @inferred forward(model, [])
        @report_call forward(model, [])
        for k in fieldnames(RockphyElastic)
            @test all(isapprox.(getfield(out_, k), getfield(outs[m], k), rtol=1e-2))
        end
    end

end
