@testitem "rock-physics tests" tags = [:rp] begin
    using MT
    methods_list = [
        SEO3,
        UHO2014, # mineral
        Ni2011,
        Sifre2014,
    ] # melt

    T = 1273.0;
    Ch2o_ol = 2e4;
    Ch2o_m = 2e4;
    Cco2_m = 200.0;

    inps = (;
        SEO3=[T],
        UHO2014=[T, Ch2o_ol],

        # melt
        Ni2011=[T, Ch2o_m],
        Sifre2014=[T, Ch2o_m, Cco2_m]
    )

    outs = (;
        SEO3=log10(8.7663e-05),
        UHO2014=log10(18.7398),

        # melt
        Ni2011=log10(0.0044),
        Sifre2014=log10(0.9710)
    )

    @testset "$m" for m in methods_list
        model = m(inps[Symbol("$m")]...)
        out_ = forward(model, [])
        @test round(out_; digits=2) ≈ round(outs[Symbol("$m")]; digits=2)
    end

    # mixing models
end