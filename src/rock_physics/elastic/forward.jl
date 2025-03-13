
function forward(m::anharmonic, p; params=default_params_anharmonic.Isaak1992)
    @unpack T_K_ref, P_Pa_ref, Gu_0_ol, dG_dT, dG_dP, ν, Gu_0_crust, dG_dT_crust, dG_dP_crust, Gu_TP, Ku_TP = params

    Gu₀, dG_dT₀, dG_dP₀ = @. calc_Gu₀(Gu_0_ol, dG_dT, dG_dP, Gu_0_crust, dG_dT_crust, dG_dP_crust) #since χ is 1., we are always using ol
    
    ΔT = @. m.T - T_K_ref # K
    ΔP = @. m.P * 1.0f9 - P_Pa_ref # Pa
    Gu_tp = @. calc_Gu(Gu₀ * 1.0f9, ΔT, ΔP, dG_dT₀, dG_dP₀, Gu_TP)
    Ku_tp = @. calc_Ku(Gu_tp, ν, Ku_TP)

    Vp = @. calc_Vp(Ku_tp, Gu_tp, m.ρ)
    Vs = @. calc_Vs(Gu_tp, m.ρ)

    return RockphyElastic(Gu_tp, Ku_tp, Vp, Vs)
end

function forward(m::anharmonic_poro, p; params=default_params_anharmonic_poro)
    @unpack m_A, m_K, ν, p_anharmonic = params

    anh_p = forward(anharmonic(m.T, m.P, m.ρ), []; params=params.p_anharmonic)

    # @show typeof(anh_p.G)

    Γ_G = @. melt_shear_moduli(m.ϕ, m_A, ν)
    Gueff = @. Γ_G * anh_p.G
    # Gueff, Γ_G = ntuple(i -> getindex.(G_s, i), length(first(G_s)))
    Γ_K = @. melt_bulk_moduli(m.ϕ, m_A, ν)
    K_sk = @. Γ_K * anh_p.K

    nr = @. (1 - K_sk / anh_p.K)^2
    dr = @. 1 - m.ϕ - K_sk / anh_p.K + m.ϕ * anh_p.K / m_K

    Kueff = @. K_sk + (nr / (dr + 1.0f-10)) * anh_p.K
    μ_eff = @. anh_p.G * Γ_G
    K_nr = @. (1 - Γ_K)^2
    K_dr = @. 1 - m.ϕ - Γ_K + m.ϕ * anh_p.K / m_K
    K_eff = @. anh_p.K * (Γ_K + K_nr / (K_dr + 1.0f-10))

    Vp = @. calc_Vp(K_eff, μ_eff, m.ρ)
    Vs = @. calc_Vs(μ_eff, m.ρ)

    # Vp_Vs = @. Vp_Vs_calc(m.ϕ, anh_p.G, anh_p.K, Γ_G, Γ_K, m.ρ, m_K) # slightly different here

    # Vp, Vs = ntuple(i -> getindex.(Vp_Vs, i), length(first(Vp_Vs)))

    return RockphyElastic(Gueff, Kueff, Vp, Vs)
end

function forward(m::SLB2005, p; params=(;))
    dV_P = @. 0.0380f0 * m.P
    dV_T = @. -0.000378f0 * (m.T - 300)

    Vs = @. 4.77 + dV_P + dV_T

    return RockphyElastic(0.0f0, 0.0f0, 0.0f0, Vs * 1.0f3)
end

function forward(::Type{M}) where {M <: AbstractElasticModel}
    return RockphyElastic
end

default_params(::Val{anharmonic}) = default_params_anharmonic.Isaak1992
default_params(::Val{anharmonic_poro}) = default_params_anharmonic_poro
default_params(::Val{SLB2005}) = (;)
