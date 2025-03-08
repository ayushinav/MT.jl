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

const params_Jones2012 = (S=10.0f0^(3.05f0), r=0.86, H=0.91, a=0.09)

const params_Poe2010 = (S_H100=10.0f0^2.59f0, H_H100=1.26f0, a_H100=1.18f0, r_H100=1.0f0,
    S_H010=10.0f0^3.46f0, H_H010=1.5f0, a_H010=1.43f0, r_H010=1.0f0,
    S_H001=10.0f0^1.02f0, H_H001=0.812f0, a_H001=0.70f0, r_H001=1.0f0, S_A100=334.0f0,
    H_A100=1.46f0, S_A010=13.8f0, H_A010=1.12f0, S_A001=99.0f0, H_A001=1.29f0)

const params_Wang2006 = (
    S_H=10.0f0^3.0f0, H_H=87.0f0, a_H=0.0f0, r_H=0.62f0, S_A=10.0f0^2.4f0, H_A=154.0f0)

const params_Yoshino2009 = (S_i=10.0f0^4.73f0, H_i=2.31f0, S_h=10.0f0^2.98f0, H_h=1.71f0,
    S_p=10.0f0^1.90f0, H_p=0.92f0, a_p=0.16f0, r_p=1.0f0)

# melts

const params_Ni2011 = (T_corr=1146.8f0, D=0.006)

const params_Sifre2014 = (
    D_p=0.007f0, D_o=0.002f0, den_p=3.3f0, den_h2o=1.4f0, den_carb=2.4f0, den_basalt=2.8f0,

    # H2O melt
    a_h2o=88774.0f0, b_h2o=0.3880f0, c_h2o=73029.0f0, d_h2o=4.54f-5, e_h2o=5.5607f0,

    # C2O melt
    a_c2o=789166.0f0, b_c2o=0.1808f0, c_c2o=32820.0f0, d_c2o=5.50f-5, e_c2o=5.7956f0)

const params_Gaillard2008 = (S=3440.0f0, H=31.9f0)
