# minerals

const params_SEO3 = (S_bfe=5.06f24, H_bfe=0.357f0, S_bmg=4.58f26, H_bmg=0.752f0,
    S_ufe=12.2f-6, H_ufe=1.05f0, S_umg=2.72f-6, H_umg=1.09f0)

const params_UHO2014 = (
    # Ionic Vacancy
    H_v=239.0f0, S_v=10.0f0^5.07f0,

    # Polaron Hopping
    H_p=144.0f0, S_p=10.0f0^2.34f0,

    # Proton
    H_h=89.0f0, S_h=10.0f0^-1.37f0, a_h=1.79f0, r_h=1.0f0)

# melts

const params_Ni2011 = (T_corr=1146.8f0, D=0.006)

const params_Sifre2014 = (
    D_p=0.007, D_o=0.002, den_p=3.3, den_h2o=1.4, den_carb=2.4, den_basalt=2.8,

    # H2O melt
    a_h2o=88774, b_h2o=0.3880, c_h2o=73029, d_h2o=4.54e-5, e_h2o=5.5607,

    # C2O melt
    a_c2o=789166, b_c2o=0.1808, c_c2o=32820, d_c2o=5.50e-5, e_c2o=5.7956)
