function forward(m::HZK2011, p; params = default_params_HZK2011)
    @unpack mechs, p_dep_calc, melt_enhancement = params

    P = @. p_dep_calc * m.P
    ϵ_rate = zero(m.T .+ m.P .+ m.dg .+ m.σ .+ m.ϕ) ## TODO
    x_ϕ_c_vec = get_melt_settings_for_x_ϕ_c(Val{melt_enhancement}())

    for mech in keys(mechs)
        sr = sr_flow_law_calculation(m.T, P* 1f9, m.σ, m.dg, m.ϕ, 0, getfield(x_ϕ_c_vec, mech), [getfield(mechs, mech)])
        @. ϵ_rate += sr
    end

    η = @. m.σ * 1f3 * 1f6 / ϵ_rate

    return RockphyViscous(ϵ_rate, η)

end

function forward(m::HK2003, p; params = default_params_HK2003)
    @unpack mechs, ch2o_o, p_dep_calc, melt_enhancement = params

    P = @. p_dep_calc * m.P

    fH2O = @. calc_fH2O(m.Ch2o, ch2o_o, P, m.T)

    ϵ_rate = zero(m.T .+ m.P.+ m.dg.+ m.σ.+ m.ϕ .+ m.Ch2o)

    # x_phi_c = Int(melt_enhancement)
    x_ϕ_c_vec = get_melt_settings_for_x_ϕ_c(Val{melt_enhancement}())
    ps_ = broadcast((x,y) ->HK2003_mech(x, y, mechs, :diff), m.T, fH2O)
    # @show ps_
    # @show Base.summarysize(mechs)
    pvec = zeros(Float32, 8);

    for mech in keys(mechs)
        ps_ = broadcast((x,y) ->HK2003_mech(x, y, mechs, mech), m.T, fH2O)
        # @show ps_[1]

        # sr = broadcast()
        # sr = broadcast((T, P, σ, d, ϕ, fH2O) -> sr_flow_law_calculation2(T, P * 1f9, σ, d, ϕ, fH2O, getfield(x_ϕ_c_vec, mech), mechs, mech, pvec), m.T, P, m.σ, m.dg, m.ϕ, fH2O)
        sr = sr_flow_law_calculation(m.T, P* 1f9, m.σ, m.dg, m.ϕ, fH2O, getfield(x_ϕ_c_vec, mech), ps_)
        # sr = sr_flow_law_calculation2(m.T, P* 1f9, m.σ, m.dg, m.ϕ, fH2O, getfield(x_ϕ_c_vec, mech), ps_)
        @. ϵ_rate += sr
    end

    η = @. m.σ * 1f3 * 1f6 / ϵ_rate

    return RockphyViscous(ϵ_rate, η)

end

function forward(m::xfit_premelt, p; params = default_params_xfit_premelt)
    @unpack α, T_η, γ, B, Tr, Pr, η_r, H, V, M, dg_r = params

    Tprime = @. m.T/m.T_solidus
    A_n = @. calc_An(Tprime, m.ϕ, α, T_η, γ, B)

    # η_meltfree = @. calc_η_meltfree(m.P, m.dg, m.T, params)

    η_meltfree = @. η_r * (m.dg/dg_r) ^ M * exp(V/(MT.gas_R * 1f3) * (m.P/m.T - Pr/ Tr)*1f9 + H/(MT.gas_R * 1f3) * (1/m.T - 1/Tr)) 

    η = @. A_n * η_meltfree

    return RockphyViscous(zero(η), η)
end