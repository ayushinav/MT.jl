using Pkg
Pkg.activate(".")

using MT
using UnPack
using SpecialFunctions
# ================= ELASTIC =================

# anharmonic

# VBR.in.SV.T_K temperature in degrees K
# VBR.in.SV.GPa pressure in GPa
# VBR.in.SV.rho density in kg m-3

params_anharmonic = (
    Isaak1992 = (
        T_K_ref = 300,
        P_Pa_ref = 1f5,
        Gu_0_ol = 81,
        dG_dT = -13.6f6,
        dG_dP = 1.8,
        ν = 0.25,
        Gu_0_crust = 40,
        dG_dT_crust = -3f6,
        dG_dP_crust = 3.5
    ),

    Cammarono2003 = (
        T_K_ref = 300,
        P_Pa_ref = 1f5,
        Gu_0_ol = 81,
        dG_dT = -14f6,
        dG_dP = 1.4,
        ν = 0.25,
        Gu_0_crust = 40,
        dG_dT__crust = -3f6,
        dG_dP_crust = 3.5
    )
)

# el_calc_Gu_0

function el_calc_Gu_0(; χ = 1f0, params_anharmonic = params_anharmonic.Isaak1992)
    
    @unpack T_K_ref, P_Pa_ref, Gu_0_ol, dG_dT, dG_dP, ν, Gu_0_crust, dG_dT_crust, dG_dP_crust = params_anharmonic

    Gu_new = (Gu_0_ol * χ + (1 - χ) * Gu_0_crust) * 1f9
    dG_dT_new = dG_dT * χ + (1 - χ) * dG_dT_crust
    dG_dP_new = dG_dP * χ + (1 - χ) * dG_dP_crust

    return (Gu_new, dG_dT_new, dG_dP_new)

end

el_calc_Gu_0()

calc_Ku(Gu, nu) = 2/3 * Gu * (1+nu)/(1-2*nu)
calc_Gu(Gu_0, dT, dP, dG_dT, dG_dP) = Gu_0 + dT*dG_dT + dP*dG_dP

function el_VpVs_unrelaxed(K, G, rho)
    Vp = sqrt((K + 4/3 * G)/rho)
    Vs = sqrt(G/rho)
    return (Vp, Vs)
end

function el_ModUnrlx_dTdP_f(T_K, P, rho, Gu_new, dG_dT_new, dG_dP_new; params_anharmonic = params_anharmonic.Isaak1992)
    @unpack T_K_ref, P_Pa_ref, Gu_0_ol, dG_dT, dG_dP, ν, Gu_0_crust, dG_dT_crust, dG_dP_crust = params_anharmonic

    dG_dT0 = dG_dT_new
    dG_dP0 = dG_dP_new
    Gu0 = Gu_new

    nu =ν

    # if default values
    # Gu_TP
    # Ku_TP
    # Gu0 = Gu_TP

    # else

    dT = T_K - T_K_ref
    dP = P* 1f9 - P_Pa_ref

    Gu_TP = calc_Gu(Gu0, dT, dP, dG_dT, dG_dP)
    Ku_TP = calc_Ku(Gu_TP, nu)

    Vp, Vs = el_VpVs_unrelaxed(Ku_TP, Gu_TP, rho)

    return (
        Gu = Gu_TP,
        Ku = Ku_TP,
        Vp = Vp,
        Vs = Vs
    )

end

el_calc_Gu_0()
el_ModUnrlx_dTdP_f(1273, 2, 3300, el_calc_Gu_0()...)

# ================= VISCOUS =================

# HZK2011
# VBR.in.SV.T_K % temperature [K]
# VBR.in.SV.P_GPa % pressure [GPa]
# VBR.in.SV.dg_um % grain size [um]
# VBR.in.SV.sig_MPa % differential stress [MPa]
# VBR.in.SV.phi % melt fraction / porosity


params_HZK2011 = (
    diff = (
        A = 10f0 ^76f-1,  # Preexponential for coble diffusion creep
        Q = 375f3,   # Activation energy for coble diffusion creep
        V = 10f-6,   # Activation volume for coble diffusion creep
        p = 3,       # Grain size exponent
        alf = 25,    # Melt factor
        r = 0,       # Water fugacity exponent
        n = 1,       # Stress exponent
        phi_c = 1f-5,
        x_phi_c = 5
    ),

    disl = (
        A = 1.1f5,   # Preexponential
        Q = 530f3,   # Activation energy
        V = 15f-6,   # Activation volume
        n = 35f-1,     # Stress exponent
        p = 0,       # Grain size exponent
        alf = 30,    # Melt factor
        r = 0,       # Water fugacity exponent
        phi_c = 1f-5,
        x_phi_c = 1
    ),

    gbs = (
        A = 10f0 ^4.8f0,  # Preexponential for GBS-disl creep
        Q = 445f3,   # Activation energy for GBS-disl creep
        V = 15f-6,   # Activation volume
        p = 73f-2,    # Grain size exponent
        n = 29f-1,     # Stress exponent
        alf = 35,    # Melt factor
        r = 0,       # Water fugacity exponent
        phi_c = 1f-5,
        x_phi_c = 25f-1
    )
)

function melt_ehancement(phi,alpha,x_phi_c,phi_c)
    a = log(x_phi_c)
    ratefac = inv(phi_c)
    step = a*erf(phi.*ratefac)
    slope = alpha*phi
    ln_SR_phi_enh = slope + step
    SR_phi_enh = exp( ln_SR_phi_enh )

    return SR_phi_enh
end

function sr_flow_law_calculation(T, P, sig, d, phi, fH2O, params)
    @unpack A, Q, V, p, n, alf, r, phi_c, x_phi_c = params

    sr = inv(x_phi_c) * A * (sig^ n) * (d^ (-p)) * exp(- (Q + P * V)/(MT.gas_R * T * 1f3)) * (fH2O^r)
    enhance = melt_ehancement(phi, alf, x_phi_c, phi_c)

    return sr * enhance
end


function HZK2011(T_K, P_GPa, dg_um, sig_MPa, phi, mechs = [:diff, :disl, :gbs]; melt_ehancement = false, )
    P_Pa = P_GPa * 1f9

    # if !(pressure dependence)
    # P_Pa = 0

    sr_tot = 0f0;

    for mech in mechs
        sr = sr_flow_law_calculation(T_K, P_Pa, sig_MPa, dg_um, phi, 0, getfield(params_HZK2011, mech))
        @show sr
        sr_tot += sr
    end

    η = sig_MPa * 1f6 / sr_tot

    return η, sr_tot

end


HZK2011(1273f0, 2f0, 4f0, 10f0, 1f-2)



