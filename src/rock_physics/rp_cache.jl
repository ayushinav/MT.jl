# minerals

const params_SEO3 = (
    S_bfe = 5.06f24, 
    H_bfe = 0.357f0,
    S_bmg = 4.58f26,
    H_bmg = 0.752f0,
    S_ufe = 12.2f-6,
    H_ufe = 1.05f0,
    S_umg = 2.72f-6,
    H_umg = 1.09f0,
)

const params_UHO2014 = (
    # Ionic Vacancy
    H_v = 239f0,         # kJ/mol, activation enthalpy
    S_v = 10f0 ^5.07f0,  # S/m, pre-exponential conductivity

    # Polaron Hopping
    H_p = 144f0,         # kJ/mol, activation enthalpy
    S_p = 10f0 ^2.34f0,  # S/m, pre-exponential conductivity
    
    # Proton
    H_h = 89f0,          # kJ/mol, activation enthalpy
    S_h = 10f0 ^-1.37f0, # S/m, pre-exponential conductivity
    a_h = 1.79f0,        # (kJ/mol/wt)*(ppm ^1/3)
    r_h = 1f0,           # unitless
)

# melts

const params_Ni2011 = (
    T_corr = 1146.8f0, # K
    D = 0.006 # unitless, Partition coefficient {ol/melt}
)

const params_Sifre2014 = (
    D_p = 0.007, # unitelss, D_{perid/melt}
    D_o = 0.002, # unitelss, D_{ol/melt}

    den_p = 3.3, # g/cm^3, density peridotite
    den_h2o = 1.4, #  g/cm^3, density of water
    den_carb = 2.4, #  g/cm^3, density of molten carbonates [Liu and Lange, 2003]
    den_basalt = 2.8, #  g/cm^3, density of molten basalt [Lange and Carmichael, 1990]

    # H2O melt
    a_h2o = 88774,
    b_h2o = 0.3880,
    c_h2o = 73029,
    d_h2o = 4.54e-5,
    e_h2o = 5.5607,

    # C2O melt
    a_c2o = 789166,
    b_c2o = 0.1808,
    c_c2o = 32820,
    d_c2o = 5.50e-5,
    e_c2o = 5.7956
)
