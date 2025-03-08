

## response

mutable struct RockphyElastic{T1,T2,T3,T4} #<: AbstractRockphyResponse
    G::T1
    K::T2
    Vp::T3
    Vs::T4
    # dG_dT::T5
    # dG_dP::T6
end

"""
T : Temperature (Kelvin)
P : Pressure (GPa)
ρ : density (kg/m³)
"""
mutable struct anharmonic{T1,T2,T3}
    T::T1
    P::T2
    ρ::T3
    # Gu_TP::T4
    # Ku_TP::T5
end

# anharmonic(T, P, ρ) = anharmonic(T, P, ρ, -1.0f0, -1.0f0)

function calc_Gu₀(G1, dG_dT1, dG_dP1, G2, dG_dT2, dG_dP2; χ=1.0f0)

    Gu₀ = (G1 * χ + (1 - χ) * G2)
    dG_dT₀ = dG_dT1 * χ + (1 - χ) * dG_dT2
    dG_dP₀ = dG_dP1 * χ + (1 - χ) * dG_dP2

    return (Gu₀, dG_dT₀, dG_dP₀)
end

calc_Ku(G, ν) = 2G/3 * (1 + ν)/ (1 - 2ν)
calc_Gu(Gu_0, ΔT, ΔP, ∂G_∂T, ∂G_∂P) = Gu_0 + ΔT * ∂G_∂T + ΔP * ∂G_∂P

calc_Vp(K, G, ρ) = sqrt((K + 4G / 3) / ρ)
calc_Vs(G, ρ) = sqrt(G / ρ)

function forward(m::anharmonic; params=params_anharmonic.Isaak1992)
    @unpack T_K_ref, P_Pa_ref, Gu_0_ol, dG_dT, dG_dP, ν, Gu_0_crust, dG_dT_crust, dG_dP_crust, Gu_TP, Ku_TP = params

    Gu₀, dG_dT₀, dG_dP₀ = calc_Gu₀(Gu_0_ol, dG_dT, dG_dP, Gu_0_crust, dG_dT_crust, dG_dP_crust) #since χ is 1., we are always using ol

    # Gu₀ to be scaled by 1f9

    if Gu_TP < 0 && Ku_TP < 0
        ΔT = m.T - T_K_ref # K
        ΔP = m.P * 1.0f9 - P_Pa_ref # Pa
        Gu_TP = calc_Gu(Gu₀ * 1.0f9, ΔT, ΔP, dG_dT₀, dG_dP₀)
        Ku_TP = calc_Ku(Gu_TP, ν)
    elseif Gu_TP < 0 && Ku_TP > 0
        @warn "Unrelaxed bulk modulus provided without unrelaxed shear modulus"
    elseif Gu_TP > 0 && Ku_TP < 0
        @warn "Unrelaxed shear modulus provided without unrelaxed bulk modulus"
    end

    Vp = calc_Vp(Ku_TP, Gu_TP, m.ρ)
    Vs = calc_Vs(Gu_TP, m.ρ)

    return RockphyElastic(Gu_TP, Ku_TP, Vp, Vs)

end

mutable struct anharmonic_poro{T1,T2,T3,T4}
    T::T1
    P::T2
    ρ::T3
    ϕ::T4
    # Gu_TP::T5
    # Ku_TP::T6
end

# anharmonic_poro(T, P, ρ, ϕ) = anharmonic_poro(T, P, ρ, ϕ, -1.0f0, -1.0f0)



function melt_shear_moduli(μ, ϕ, A, ν)
    ψ = 1 - A * sqrt(ϕ)

    b = reshape(Float32.([1.6122, 4.5869, -7.5395, 0.13572, 3.6086, -4.8676, 0, 0, -4.3182]), 3, 3)
    ν_vec = ν .^ [0, 1, 2]

    b_vec = b * ν_vec # :)

    n_μ = b_vec[1] * ψ + b_vec[2] * (1 - ψ) + b_vec[3] * ψ * (1 - ψ)^2

    μ_sk_prime = (1 - (1 - ψ)^n_μ)
    Γ_G = (1 - ϕ) * μ_sk_prime
    μ_sk = Γ_G * μ

    return μ_sk, Γ_G
end

function melt_bulk_moduli(k, ϕ, A, Km, ν)
    ψ = 1f0 - A * sqrt(ϕ)

    a = reshape(Float32.([1.8625, 4.5001, -5.6512, 0.52594, -6.1551, 6.9159, -4.8397, -4.3634, 29.595, 0.0, 0.0, -58.96]), 3, 4)
    ν_vec = ν .^ [0, 1, 2, 3]

    a_vec = a * ν_vec # :)

    n_k = a_vec[1] * ψ + a_vec[2] * (1 - ψ) + a_vec[3] * ψ * (1 - ψ)^(1.5f0)

    k_sk_prime = (1 - (1 - ψ)^n_k)
    Γ_K = (1 - ϕ) * k_sk_prime
    K_sk = Γ_K * k

    nr = (1 - K_sk/k)^2
    dr = 1 - ϕ - K_sk/ k + ϕ*k/ Km

    KB_eff = K_sk + (nr / (dr + 1.0f-10)) * k

    return KB_eff, Γ_K
end

function Vp_Vs_calc(ϕ, G, K, Γ_G, Γ_K, ρ, Km)
    μ_eff = G * Γ_G

    K_nr = (1- Γ_K)^2
    K_dr = 1 - ϕ - Γ_K + ϕ *K/Km

    K_eff = K * (Γ_K + K_nr/(K_dr + 1f-10))

    Vp = calc_Vp(K_eff, μ_eff, ρ)
    Vs = calc_Vs(μ_eff, ρ)

    return Vp, Vs    
end

function forward(m::anharmonic_poro; params=params_anharmonic_poro)
    @unpack m_A, m_K, ν, p_anharmonic = params

    anh_p = forward(
        anharmonic(m.T, m.P, m.ρ); params=params.p_anharmonic
    )

    Gueff, Γ_G = melt_shear_moduli(anh_p.G, m.ϕ, m_A, ν)
    Kueff, Γ_K = melt_bulk_moduli(anh_p.K, m.ϕ, m_A, m_K, ν)

    Vp, Vs = Vp_Vs_calc(m.ϕ, anh_p.G, anh_p.K, Γ_G, Γ_K, m.ρ, m_K) # slightly different here

    return RockphyElastic(Gueff, Kueff, Vp, Vs)
end


mutable struct SLB2005_{T1,T2}
    T::T1
    P::T2
end


function forward(m::SLB2005_; params = [])
    dV_P = 0.0380f0 * m.P
    dV_T = -0.000378f0 * (m.T- 300)

    Vs = 4.77 + dV_P + dV_T

    return RockphyElastic(0f0, 0f0, 0f0, Vs * 1f3)
end

m_anharmonic = anharmonic(1273f0, 0.2f0, 3300f0)

forward(m_anharmonic)

m_anharmonic_poro = anharmonic_poro(1273f0, 0.2f0, 3300f0, 0.01f0)

forward(m_anharmonic_poro)

m_SLB2005 = SLB2005_(1273f0, 0.2f0)

forward(m_SLB2005)

