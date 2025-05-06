default_params_Hirschmann2000 = (
    A1 = 1108.08f0, # [C]
    A2 = 139.44f0, # [C/GPa]
    A3 = -5.904f0 # [C/GPa^2]
)

"""
returns dry solidus according to Hirschmann2000

## Arguments
 - `P` : GPa
"""
function solidus_Hirschmann2000(ps_nt, params = default_params_Hirschmann2000)
    @unpack A1, A2, A3 = params
    @unpack P = ps_nt
    T_solidus = @. A1 + A2 * P + A3 * P*P
    return (; T_solidus)
end

"""
returns pressure-dependent dry solidus according to Katz2003

## Arguments
 - `P` : GPa
"""
function solidus_Katz2003(ps_nt, params = default_params_Katz2003)
    @unpack A1, A2, A3 = params
    @unpack P = ps_nt
    T_solidus = @. A1 + A2 * P + A3 *P*P
    return (; T_solidus)
end

"""
updates solidus with the depressed one according to Katz2003

## Arguments
 - `T` : solidus temperature (K)
 - `Ch2o_m` : melt fraction
"""
function supress_solidus_Katz2003(ps_nt)
    @unpack Ch2o_m, P, T_solidus = ps_nt
    γ = 0.75f0
    K = 43f0

    Ch2o_m_sat = @. 12 * P^(0.6f0) + P
    # C_w = (Ch2o_m <= Ch2o_m_sat) ? (Ch2o_m) : Ch2o_m # worth investigating
    Ch2o_m_ = @. Ch2o_m * 1f-4
    dT = @. K * (Ch2o_m_)^γ
    T_solidus = @. T_solidus - dT
    return (; T_solidus)
end

function Dasgupta2007_core(Cco2_m_)

    if Cco2_m_ <= 25
        dT = 27.04f0 * Cco2_m_ +  1490.75f0 * log((100f0-1.18f0 * Cco2_m_)/100f0);
    elseif Cco2_m_ > 25 && Cco2_m_ <= 37
        dTmax = 27.04f0 * 25f0 +  1490.75f0 * log((100f0-1.18f0 * 25f0)/100f0);
        dT = dTmax +  (160f0 - dTmax)/(37f0 - 25f0) * (Cco2_m_-25f0);
    elseif Cco2_m_ > 37
        dTmax = 27.04f0 * 25f0 +  1490.75f0 * log((100f0-1.18f0 * 25f0)/100f0);
        dTmax = dTmax +  (160f0 - dTmax);
        dT = dTmax + 150f0;
    end

    return dT
end

function supress_solidus_Dasgupta2007(ps_nt)
    @unpack Cco2_m, T_solidus = ps_nt

    Cco2_m_ = @. Cco2_m * 1f-4
    dT = @. Dasgupta2007_core(Cco2_m_)

    T_solidus = @. T_solidus - dT
    return (; T_solidus)
end

function supress_solidus_Blatter2022(ps_nt)
    @unpack Cco2_m, Ch2o_m, T_solidus, D = ps_nt

    # D = 0.005
    Ch2o_m_ = @. Ch2o_m/D * 1f-4

    M = 59f0
    ΔS = 0.4

    X_OH = @. 2M*Ch2o_m_/100 / (18.02 + Ch2o_m_/100 * (2M - 18.02))
    # X_OH = 2M*Ch2o_m_ / (18.02 + Ch2o_m_ * (2M - 18.02))
    # @show X_OH, Ch2o_m_

    T_wet = @. T_solidus* inv(1 - MT.gas_R * 1f3/(M * ΔS) * log(1 - X_OH))

    Cco2_m_ = @.Cco2_m * 1f-4

    #=
    2 GPa: a = 19.21; b = 1491.37; c = 0.86 (this study) 
    3 GPa: a = 27.04; b = 1490.75; c = 1.18 (ref. 5) 
    4 GPa: a = 31.90; b = 1469.92; c = 1.31 (this study) 
    5 GPa: a = -5.01; b = 1514.84; c = -1.23 (this study) 
    =#

    a = 19.21f0
    b = 1491.37f0
    c = 0.86f0

    ΔT_co2 = @. a * Cco2_m_ + b * log(1 - c * Cco2_m_ * 1f-2)
    T_solidus = @. T_wet - ΔT_co2
    return (; T_solidus)
end

function get_Cco2_m(ps_nt)
    @unpack ϕ, Cco2 = ps_nt
    Cco2_m = @. Cco2 * inv(ϕ)
    # Cco2_ol = @. Cco2 * inv(1- ϕ)
    # Cco2_ol = @. Cco2 - Cco2_m
    return (; Cco2_m)
end

function get_Ch2o_m(ps_nt)
    @unpack ϕ, Ch2o, D = ps_nt
    Ch2o_m = @. Ch2o * inv(D + ϕ* (1-D))
    # Ch2o_ol = @. Ch2o * inv((1-D) + (1-ϕ)* (D))
    # Ch2o_ol = @. Ch2o - Ch2o_m
    return (; Ch2o_m)
end


function get_melt_fraction(ps_nt)
    @unpack Cco2, Ch2o, T_solidus, P, D = ps_nt

    ϕ = broadcast(
        get_melt_fraction_core, Cco2, Ch2o, T_solidus, P, D
    )

    return (; ϕ)
end


function get_melt_fraction_core(Ch2o, Cco2, T_solidus, P, D)
    function f(u, p)

        Ch2o_m = get_Ch2o_m((;ϕ = u, p.Ch2o, p.D)).Ch2o_m
        Cco2_m = get_Cco2_m((;ϕ = u, p.Cco2)).Cco2_m

        # T_new = supress_solidus_Blatter2022((;Ch2o_m, Cco2_m, p.T_solidus)).T_solidus
        T_new1 = supress_solidus_Dasgupta2007((;Ch2o_m, Cco2_m, p.T_solidus)).T_solidus
        T_new2 = supress_solidus_Katz2003((;Ch2o_m, Cco2_m, p.T_solidus, p.P)).T_solidus
        ΔT = 2p.T_solidus - T_new1 - T_new2
        dTdF = -40*p.P + 450

        # return u - ΔT/dTdF
        return u * dTdF/ΔT - 1
    end

    prob_init = IntervalNonlinearProblem(f, (1f-15, 1f0), (; Ch2o, Cco2, T_solidus, P, D))
    sol = solve(prob_init)
    return sol.u
end


function get_melt_fraction(ps_nt)
    @unpack Cco2, Ch2o, T_solidus, P, D = ps_nt

    ϕ = broadcast(
        get_melt_fraction_core, Cco2, Ch2o, T_solidus, P, D
    )

    return (; ϕ)
end