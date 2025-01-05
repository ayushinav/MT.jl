@testitem "rock-physics tests" tags = [:rp] begin
    using MT
    methods_list = [SEO3, UHO2014, # mineral
        Ni2011, Sifre2014]; # melt

    inps = (; zip(
        Symbol.(methods_list), [
            [1273.0], # SEO3
            [1273.0, 4000.0], # UHO2014
            [1273.0, 400.0], # Ni2011
            [1273.0, 400.0, 400.0], # Sifre2014
        ]
    )...)

    outs =
        (; zip(
            Symbol.(methods_list),
            log10.([9.0720e-04, 0.0017, 0.4185, 0.7552])
        )...)


    @testset "$m" for m in methods_list
        # @show m, typeof(m)
        # @show inps[Symbol("$m")]
        model = m(inps[Symbol("$m")]...)
        out_ = forward(model, [])
        @test out_ == outs[Symbol("$m")]
    end

    # UHO failing

    # mixing models
end
