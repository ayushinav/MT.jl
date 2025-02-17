@testitem "rock-physics tests" tags = [:rp] begin
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
        Gaillard2008
    ]

    T = 1273.0
    Ch2o_ol = 2e4
    Ch2o_m = 2e4
    Cco2_m = 200.0

    inps = (;
        SEO3 = [T],
        UHO2014 = [T, Ch2o_ol],
        Jones2012 = [T, Ch2o_ol],
        Yoshino2009 = [T, Ch2o_ol],
        Wang2006 = [T, Ch2o_ol],
        Poe2010 = [T, Ch2o_ol],

        # melt
        Ni2011=[T, Ch2o_m],
        Sifre2014=[T, Ch2o_m, Cco2_m],
        Gaillard2008=[T]
    )

    outs = (;
        SEO3 = log10(8.7663e-05),
        UHO2014 = log10(18.7398),
        Jones2012 = log10(1.4290),
        Yoshino2009 = log10(0.2276),
        Wang2006 = log10(0.4138),
        Poe2010 = log10(2.8011e3),

        # melt
        Ni2011=log10(0.0044),
        Sifre2014=log10(0.9710),
        Gaillard2008=log10(168.8759)
    )

    @testset "$(methods_list[i])" for i in eachindex(methods_list)
        m = keys(outs)[i]
        model = methods_list[i](inps[m]...)
        out_ = forward(model, [])
        @test round(out_; digits = 2) ≈ round(outs[m]; digits = 2)
    end

    # mixing models
    mixing_list = [HS1962_plus(), HS1962_minus(), MAL(0.2)]
    mixing_outs = log10.([3.869f-4, 1.1502f-4, 2.8f-3])

    @testset "$(mixing_list[i])" for i in eachindex(mixing_list)
        model = construct_mixing_models([1000.0 + 273.0, 2e4],
            [:T, :Ch2o_m],
            [0.1],
            [SEO3, Ni2011],
            [mixing_list[i]])

        out_ = forward(model, [])
        @test round(first(out_.σ); digits=2) ≈ round(mixing_outs[i]; digits=2)
    end
end
