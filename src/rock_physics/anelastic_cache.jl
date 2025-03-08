params_andrade_psp = (
    n = inv(3f0),
    β = 2f-2,
    τ_MR = 10f0 ^ (5.3f0),
    E = 303f3,
    G_UR = 62.2f0, # GPa,
    TR = 1173,
    PR = 2f-1,
    dR = 3.1f0,
    Vstar = 1f-5,
    M = 1,
    melt_alpha = 25,
    ϕ_c = 1f-5,
    # x_ϕ_c = 5,
    elastic_type = anharmonic,
    params_elastic = params_anharmonic.Isaak1992,
    melt_enhancement = false
)


params_JF10 = (
    bg_only = (
        dR = 13.4,  # ref grain size in microns
        G_UR = 62.5,  # GPa, unrel. G, reference val.
        E = 303000,  # J/mol
        m_a = 1.19,  # grain size exponent for tau_i, i in (L,H,P)
        alf = 0.257,  # high temp background tau exponent
        DeltaB = 1.13,  # relaxation strength
        Tau_LR = 1e-3,  # Relaxation time lower limit reference
        Tau_HR = 1e7,  # Relaxation time higher limit reference
        Tau_MR = 10^6.95,  # Reference Maxwell relaxation time
        DeltaP = 0,  # no peak, set to 0
        sig = 0,  # no peak, set to 0
        Tau_PR = 0,  # no peak, set to 0
        TR = 1173,
        PR = 0.2,
        Vstar = 1f-5,
        m_v = 3, 
        melt_alpha = 25,
        ϕ_c = 1f-5,
        # x_ϕ_c = 5,
    ),
    bg_peak = (
        DeltaP = 0.057,  # relaxation strength of peak
        sig = 4,  # sigma, peak breadth
        Tau_PR = 10^-3.4,  # center maxwell time
        dR = 13.4,  # ref grain size in microns
        G_UR = 66.5,  # GPa, unrel. G, reference val.
        E = 360000,  # J/mol
        m_a = 1.31,  # grain size exponent for tau_i, i in (L,H,P)
        alf = 0.274,  # high temp background tau exponent
        DeltaB = 1.13,  # relaxation strength of background
        Tau_LR = 1e-3,  # Relaxation time lower limit reference
        Tau_HR = 1e7,  # Relaxation time higher limit reference
        Tau_MR = 10^7.48,  # Reference Maxwell relaxation time
        TR = 1173,
        PR = 0.2,
        Vstar = 1f-5,
        m_v = 3,
        melt_alpha = 25,
        ϕ_c = 1f-5,
        # x_ϕ_c = 5,
    ),
    s6585_bg_only = (
        dR = 3.1,  # ref grain size in microns
        G_UR = 62.0,  # GPa, unrel. G, reference val.
        E = 303000,  # J/mol
        m_a = 1.19,  # grain size exponent for tau_i, i in (L,H,P)
        alf = 0.33,  # high temp background tau exponent
        DeltaB = 1.4,  # relaxation strength
        Tau_LR = 1e-2,  # Relaxation time lower limit reference
        Tau_HR = 1e6,  # Relaxation time higher limit reference
        Tau_MR = 10^5.2,  # Reference Maxwell relaxation time
        DeltaP = 0,  # no peak, set to 0
        sig = 0,  # no peak, set to 0
        Tau_PR = 0,  # no peak, set to 0
        TR = 1173,
        PR = 0.2,
        Vstar = 1f-5,
        m_v = 3,
        melt_alpha = 25,
        ϕ_c = 1f-5,
        # x_ϕ_c = 5,
    ),
    s6585_bg_peak = (
        DeltaP = 0.07,  # relaxation strength of peak
        sig = 4,  # sigma, peak breadth
        Tau_PR = 10^-2.9,  # center maxwell time
        dR = 3.1,  # ref grain size in microns
        G_UR = 66.5,  # GPa, unrel. G, reference val.
        E = 327000,  # J/mol
        m_a = 1.19,  # grain size exponent for tau_i, i in (L,H,P)
        alf = 0.33,  # high temp background tau exponent
        DeltaB = 1.4,  # relaxation strength
        Tau_LR = 1e-2,  # Relaxation time lower limit reference
        Tau_HR = 1e6,  # Relaxation time higher limit reference
        Tau_MR = 10^5.4,  # Reference Maxwell relaxation time
        TR = 1173,
        PR = 0.2,
        Vstar = 1f-5,
        m_v = 3,
        melt_alpha = 25,
        ϕ_c = 1f-5,
        # x_ϕ_c = 5,
    )
)

params_eburgers_psp = (
    integration_params = (
        nτ = 3000,
        integration_method = :quadgk,
        τ_integration_points = 5000
    ),

    params_btype = params_JF10.bg_only,
    elastic_type = anharmonic,
    params_elastic = params_anharmonic.Isaak1992,
    viscous_type = xfit_premelt,
    params_viscous = params_xfit_premelt,
    melt_enhancement = false,
    JF10_visc = true

)

params_premelt_anelastic = (
    params_xfit = ( 
        α_B = 0.38,  # high temp background exponent
        A_B = 0.664,  # high temp background dissipation strength

        # Pre-melting dissipation peak settings
        τ_pp = 6e-5,  # peak center, table 4 of YT16, paragraph before eq 10
        A_p_fac_1 = 0.01,
        A_p_fac_2 = 0.4,
        A_p_fac_3 = 0.03,
        σ_p_fac_1 = 4,
        σ_p_fac_2 = 37.5,
        σ_p_fac_3 = 7,
        A_p_Tn_pts = [0.91, 0.96, 1],  # Tn cutoff points
        σ_p_Tn_pts = [0.92, 1],  # Tn cutoff points

        # Melt effects
        include_direct_melt_effect = false,  # set to 1 to include YT2024 melt effect
        β = 1.38,  # determined in YT2024, named Beta_P in YT2024 eq 5
        β_B = 6.94,  # YT2024 only
        poro_Λ = 4.0 , # Table 6 YT2024
    ),

    elastic_type = anharmonic,
    elastic_params = params_anharmonic.Isaak1992,

    viscous_type = xfit_premelt,
    viscous_params = params_xfit_premelt
)

params_xfit_mxw = (
    fit1 = (
        β2 = 1853,
        τ_cutoff = 1f-11,
    ),

    fit2 = (
        β2_fit2 = 8.476f0,
        τ_cutoff_fit2 = 5f-6
    ),

    α2 = 0.5f0,
    β1 = 0.32f0,
    α_a = 0.39f0,
    α_b = 0.28f0,
    α_c = 2.6f0,
    α_τn = 1f-1,
    melt_alpha = 25,
    ϕ_c = 1f-5,
    fit = :fit1,
    scaling_method = 1,

    elastic_type = anharmonic,
    elastic_params = params_anharmonic.Isaak1992,

    viscous_type = xfit_premelt,
    viscous_params = params_xfit_premelt
)
