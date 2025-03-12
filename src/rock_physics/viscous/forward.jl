# HZK2011

function get_melt_enhancement(phi,α,x_ϕ_c,ϕ_c)
    a = log(x_ϕ_c)
    ratefac = inv(ϕ_c)
    step = a*erf(phi.*ratefac)
    slope = α*phi
    ln_SR_phi_enh = slope + step
    SR_phi_enh = exp(ln_SR_phi_enh)

    return SR_phi_enh
end


function sr_flow_law_calculation(T, P, σ, d, ϕ, fH2O, params)
    @unpack A, Q, V, p, n, alf, r, ϕ_c, x_ϕ_c = params

    # (x_phi_c == 0) && (x_ϕ_c = 1f0)

    sr = inv(x_ϕ_c) * A * ((σ* 1f3)^ n) * (d^ (-p)) * exp(- (Q + P * V)/(MT.gas_R * T * 1f3)) * (fH2O^r)
    enhance = get_melt_enhancement(ϕ, alf, x_ϕ_c, ϕ_c)

    return sr * enhance
end


function forward(m::HZK2011; params = params_HZK2011)
    @unpack mechs, p_dep_calc, melt_enhancement = params

    P = p_dep_calc * m.P

    ϵ_rate = 0f0

    # x_phi_c = Int(melt_enhancement)
    x_ϕ_c_vec = get_melt_settings_for_x_ϕ_c(Val{melt_enhancement}())

    for mech in keys(mechs)
        sr = sr_flow_law_calculation(m.T, P* 1f9, m.σ, m.dg, m.ϕ, 0, (getfield(mechs, mech)..., x_ϕ_c = getfield(x_ϕ_c_vec, mech)))
        ϵ_rate += sr
    end

    η = m.σ * 1f3 * 1f6 / ϵ_rate

    return RockphyViscous(ϵ_rate, η)

end

# HK2003

function calc_fH2O(H2O_ppm, H2O_o, P, T)
    E = 40f3
    V = 10f-6
    A_o = 26

    return (H2O_ppm >= H2O_o) * H2O_ppm/A_o * exp((E + P *V)/ (MT.gas_R * 1f3 * T))
end

function HK2003_mech(T, fH2O, mechs, mech)
    if mech == :diff
        if fH2O > 0
            ps =  getfield(getfield(mechs, mech), :wet)
        else
            ps =  getfield(getfield(mechs, mech), :dry)
        end
    elseif mech == :disl
        if fH2O > 0
            ps =  getfield(getfield(mechs, mech), :wet)
        else
            ps =  getfield(getfield(mechs, mech), :dry)
        end
    elseif mech == :gbs
        if T >= (1250 + 273)
            ps =  getfield(getfield(mechs, mech), :gt1250)
        else
            ps =  getfield(getfield(mechs, mech), :lt1250)
        end
    end

    # x_ϕ_c_vec = get_melt_settings_for_x_ϕ_c(Val{melt_enhancement}())

    return ps #(ps..., x_ϕ_c = getfield(x_ϕ_c_vec, mech))
end

function forward(m::HK2003; params = params_HK2003)
    @unpack mechs, ch2o_o, p_dep_calc, melt_enhancement = params

    P = p_dep_calc * m.P

    fH2O = calc_fH2O(m.Ch2o, ch2o_o, P, m.T)

    ϵ_rate = 0f0

    # x_phi_c = Int(melt_enhancement)
    x_ϕ_c_vec = get_melt_settings_for_x_ϕ_c(Val{melt_enhancement}())

    for mech in keys(mechs)
        sr = sr_flow_law_calculation(m.T, P* 1f9, m.σ, m.dg, m.ϕ, fH2O, (HK2003_mech(m.T, fH2O, mechs, mech)..., x_ϕ_c = getfield(x_ϕ_c_vec, mech)))
        ϵ_rate += sr
    end

    η = m.σ * 1f3 * 1f6 / ϵ_rate

    return RockphyViscous(ϵ_rate, η)

end

# xfit_premelt

function calc_An(Tprime, ϕ, params) 
    @unpack α, T_η, γ, B, Tr, Pr, η_r, H, V, M, dg_r = params

    if Tprime < T_η
        return An = 1
    else 
        if (Tprime < 1)
            return An = exp(- log(γ) *(Tprime - T_η)/ (Tprime * (1 - T_η)))
        else
            return An = exp(-α * ϕ) * inv(γ * B)
        end
    end
    
end

function calc_η_meltfree(P, dg, T, params)
    @unpack α, T_η, γ, B, Tr, Pr, η_r, H, V, M, dg_r = params

    η = η_r * (dg/dg_r) ^ M * exp(V/(MT.gas_R * 1f3) * (P/T - Pr/ Tr)*1f9 + H/(MT.gas_R * 1f3) * (1/T - 1/Tr)) 
    return η
end

function forward(m::xfit_premelt; params = params_xfit_premelt)

    Tprime = m.T/m.T_solidus
    A_n = calc_An(Tprime, m.ϕ, params)

    η_meltfree = calc_η_meltfree(m.P, m.dg, m.T, params)

    η = A_n * η_meltfree

    return RockphyViscous(0f0, η)
end