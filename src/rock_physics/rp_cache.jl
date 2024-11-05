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
    H_v = 239f0; # kJ/mol, activation enthalpy
    S_v = 10f0 ^5.07f0; # S/m, pre-exponential conductivity

    # Polaron Hopping
    H_p = 144f0; # kJ/mol, activation enthalpy
    S_p = 10f0 ^2.34f0; # S/m, pre-exponential conductivity
    
    # Proton
    H_h = 89f0; # kJ/mol, activation enthalpy
    S_h = 10f0 ^-1.37f0; # S/m, pre-exponential conductivity
    a_h = 1.79f0; # (kJ/mol/wt)*(ppm ^1/3)
    r_h = 1f0; # unitless
)

# melts

const params_Ni2011 = (
    T_corr = 1146.8f0, # K
    D = 0.006 # unitless, Partition coefficient {ol/melt}
)
