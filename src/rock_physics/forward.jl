# utils

fO2(T) = 10.0f0^(-24441.9f0 * inv(T) + 13.296f0)

# minerals

function forward(m::SEO3, p)
    @unpack S_bfe, H_bfe, S_bmg, H_bmg, S_ufe, H_ufe, S_umg, H_umg = params_SEO3
    fO₂ = fO2(m.T)

    # sT 
    bfe = S_bfe * exp(-H_bfe * inv(boltz_k * m.T))
    bmg = S_bmg * exp(-H_bmg * inv(boltz_k * m.T))
    ufe = S_ufe * exp(-H_ufe * inv(boltz_k * m.T))
    umg = S_umg * exp(-H_umg * inv(boltz_k * m.T))
    concFe = bfe + 3.33f24 * exp(-0.02f0 * inv(boltz_k * m.T)) * fO₂^(1 / 6)
    concMg = bmg + 6.21f30 * exp(-1.83f0 * inv(boltz_k * m.T)) * fO₂^(1 / 6)
    σ = concFe * ufe * charge_e + 2.0f0 * concMg * umg * charge_e

    return log10(σ)
end

function forward(m::UHO2014, p)
    @unpack H_v, S_v, H_p, S_p, H_h, S_h, a_h, r_h = params_UHO2014

    σ_v = S_v * exp(-H_v * inv(gas_R * m.T))
    σ_p = S_p * exp(-H_p * inv(gas_R * m.T))

    H_wet = H_h - a_h * (m.Ch2o_ol^inv(3.0f0))
    S_wet = S_h * (m.Ch2o_ol^r_h)
    σ_h = S_wet * exp(-H_wet * inv(gas_R * m.T))

    σ = σ_v + σ_p + σ_h

    return log10(σ)
end

# melt

function forward(m::Ni2011, p)
    @unpack T_corr, D = params_Ni2011

    ls = 2.172f0 - (860.82f0 - 204.46f0 * sqrt(m.Ch2o_m / 1.0f4)) * inv(m.T - T_corr)
    σ = 10.0f0^ls

    return log10(σ)
end

function forward(m::Sifre2014, p)
    @unpack a_h2o, b_h2o, c_h2o, d_h2o, e_h2o, a_c2o, b_c2o, c_c2o, d_c2o, e_c2o = params_Sifre2014

    H_h2o = a_h2o * exp(-b_h2o * m.Ch2o_m * 1.0f-4) + c_h2o
    lS_h2o = d_h2o * H_h2o + e_h2o
    S_h2o = exp(lS_h2o)
    melt_h2o = S_h2o * exp(-H_h2o * inv(1f3 * gas_R * m.T))

    # Esig CO2 melt
    H_co2 = a_c2o * exp(-b_c2o * m.Cco2_m * 1.0f-4) + c_c2o
    lS_co2 = d_c2o * H_co2 + e_c2o
    S_co2 = exp(lS_co2)
    melt_co2 = S_co2 * exp(-H_co2 * inv(1f3 * gas_R * m.T))

    # summation of conduction mechanisms
    σ = melt_co2 + melt_h2o

    return log10(σ)
end
