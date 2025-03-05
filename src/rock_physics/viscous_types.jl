
mutable struct RockphyViscous{T1,T2,T3,T4} <: AbstractRockphyResponse
    ϵ_rate::T1
    η::T2
end

"""
T : K
P : GPa
dg : μm
σ : GPa
ϕ
"""
mutable struct HZK2011{T1, T2, t3, T4, T5}
    T::T1,
    P::T2,
    dg::T3,
    σ::T4,
    ϕ::T5
end


function melt_ehancement(phi,α,x_ϕ_c,ϕ_c)
    a = log(x_ϕ_c)
    ratefac = inv(ϕ_c)
    step = a*erf(phi.*ratefac)
    slope = α*phi
    ln_SR_phi_enh = slope + step
    SR_phi_enh = exp( ln_SR_phi_enh )

    return SR_phi_enh
end


function sr_flow_law_calculation(T, P, σ, d, ϕ, x_phi_c, fH2O, params)
    @unpack A, Q, V, p, n, alf, r, ϕ_c, x_ϕ_c = params

    (x_phi_c == 1) && x_ϕ_c = x_phi_c

    sr = inv(x_ϕ_c) * A * (σ^ n) * (d^ (-p)) * exp(- (Q + P * V)/(MT.gas_R * T * 1f3)) * (fH2O^r)
    enhance = melt_ehancement(ϕ, alf, x_ϕ_c, ϕ_c)

    return sr * enhance
end


function forward(m::HZK2011; params = params_HZK2011)
    @unpack mechs, p_dep_calc, melt_enhancement = params

    P = p_dep_calc * m.P

    ϵ_rate = 0f0

    
    x_phi_c = Int(melt_enhancement)

    # (melt_enhancement) && ()

    for mech in (mechs)
        sr = sr_flow_law_calculation(m.T, m.P* 1f9, σ, dg, ϕ, x_phi_c, 0, getfield(params_HZK2011, mech))
        ϵ_rate += sr
    end

    η = σ * 1f-3 / ϵ_rate

    return RockphyViscous(ϵ_rate, η)

end