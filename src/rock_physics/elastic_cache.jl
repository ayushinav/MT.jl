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

params_anharmonic_poro = (
    m_A = 1f6,
    m_K = 30f9,
    ν = 0.25,
    p_anharmonic = params_anharmonic.Isaak1992
)