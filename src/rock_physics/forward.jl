# utils

fO2(T) = 10.0f0^(-24441.9f0 * inv(T) + 13.296f0)

# minerals

function forward(m::SEO3, p)
    @unpack S_bfe, H_bfe, S_bmg, H_bmg, S_ufe, H_ufe, s_umg, H_umg = params_SEO3
    fO₂ = fO2(m.T)

    # sT 
    bfe = S_bfe * exp(-H_bfe * inv(boltz_k * m.T))
    bmg = S_bmg * exp(-H_bmg * inv(boltz_k * m.T))
    ufe = S_ufe * exp(-H_ufe * inv(boltz_k * m.T))
    umg = S_umg * exp(-H_umg * inv(boltz_k * m.T))
    concFe = bfe + 3.33e24 * exp(-0.02 * inv(boltz_k * m.T)) * fO₂^(1 / 6)
    concMg = bmg + 6.21e30 * exp(-1.83 * inv(boltz_k * m.T)) * fO₂^(1 / 6)
    σ = concFe * ufe * charge_e + 2.0 * concMg * umg * charge_e

    return σ
end

function forward(m::UHO2014, p)
    @unpack H_v, S_v, H_p, S_p, H_h, S_h, a_h, r_h = params_UHO2014

    σ_v = S_v * exp(-H_v * inv(gas_R * m.T))
    σ_p = S_p * exp(-H_p * inv(gas_R * m.T))

    H_wet = H_h - a_h * (Ch2o_ol^inv(3.0f0))
    S_wet = S_h * (Ch2o_ol^r)
    σ_h = S_wet * exp(-H_wet * inv(gas_R * m.T))

    σ = σ_v + σ_p + σ_h

    return σ
end

# melt

function forward(m::Ni2011, p)
    @unpack T_corr, D = params_Ni2011

    ls = 2.172 - (860.82 - 204.46 * sqrt(m.Ch2o_m / 1e4)) * inv(m.T - Tcorr)
    σ = 10.0^ls

    return σ
end
