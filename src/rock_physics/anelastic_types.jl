mutable struct RockPhyAnelastic{T1, T2, T3, T4, T5, T6}
    J1::T1
    J2::T2
    Qinv::T3
    M::T4
    V::T5
    Vave::T6
end

mutable struct andrade_psp{T1, T2, T3, T4, T5, T6, T7, T8, T9}
    T::T1
    P::T2
    dg::T3
    σ::T4
    ϕ::T5
    ρ::T6
    Ch2o::T7
    T_solidus::T8
    f::T9
end

andrade_psp(T, P, dg, σ, ϕ, ρ, f) = andrade_psp(T, P, dg, σ, ϕ, ρ, 0f0, 0f0, f)

function calc_X̃(T, d, P, ϕ, params_anelastic)
    @unpack n, β, τ_MR, E, G_UR, TR, PR, dR, Vstar, M, melt_alpha, ϕ_c , elastic_type, params_elastic, melt_enhancement = params_anelastic

    X̃ = (d/ dR)^(-M) * exp((-E/(MT.gas_R * 1f3)) * (inv(T) - inv(TR)) -Vstar/(MT.gas_R * 1f3) * (P/T - PR/TR) * 1f9)
    x_ϕ_c = getfield(get_melt_settings_for_x_ϕ_c(Val{melt_enhancement}()), :diff)
    X̃ = X̃/x_ϕ_c
    X̃ *= get_melt_enhancement(ϕ, melt_alpha, x_ϕ_c, ϕ_c)

    return X̃

end

function forward(m::andrade_psp; params = params_andrade_psp)
    @unpack n, β, τ_MR, E, G_UR, TR, PR, dR, Vstar, M, melt_alpha, ϕ_c , elastic_type, params_elastic, melt_enhancement = params_andrade_psp

    # @show [getfield(m, k) for k in fieldnames(elastic_type)]
    resp_elastic = forward(elastic_type([getfield(m, k) for k in fieldnames(elastic_type)]...), params = params_elastic)
    # resp_elastic = forward(anharmonic(m.T, m.P, m.ρ))
    @unpack G, K, Vp, Vs = resp_elastic

    Ju = inv(G)
    ω = 2π .* m.f
    # τ₀ = inv.(f)

    X̃ = calc_X̃(m.T, m.dg, m.P, m.ϕ, params)

    param1 = β * gamma(1+n)*cos(π/2*n)
    param2 = β * gamma(1+n)*sin(π/2*n)

    ωX = ω./X̃

    J1 = Ju * (1 + param1 .* inv.(ωX).^n)
    J2 = Ju * (param2 .* inv.(ωX).^n .+ inv.(τ_MR .* ωX))
    Qinv = J2 .* inv.(J1)
    Ma = sqrt.(inv.(J1.^2 + J2.^2))
    Va = sqrt.(Ma ./ m.ρ)

    Vave = sum(Va) * inv(length(m.f))

    return RockPhyAnelastic(
        J1, J2, Qinv, Ma, Va, Vave
    )

end

mutable struct eburgers_psp{T1, T2, T3, T4, T5, T6, T7, T8, T9}
    T::T1
    P::T2
    dg::T3
    σ::T4
    ϕ::T5
    ρ::T6
    Ch2o::T7
    T_solidus::T8
    f::T9
end

eburgers_psp(T, P, dg, σ, ϕ, ρ, f) = eburgers_psp(T, P, dg, σ, ϕ, ρ, 0f0, 0f0, f)

function add_melt_effects(ϕ, scale, m_α, ϕ_c, x_ϕ_c)
    return scale * x_ϕ_c / get_melt_enhancement(ϕ, m_α, x_ϕ_c, ϕ_c)
end

function get_η_diff(m, viscous_type::Val{HZK2011}, params_viscous)
    @unpack mechs, p_dep_calc, melt_enhancement = params_viscous

    P = p_dep_calc * m.P
    x_ϕ_c = get_melt_settings_for_x_ϕ_c(Val{melt_enhancement}())
    fH2O = 0
    ϵ_rate_diff = sr_flow_law_calculation(m.T, P* 1f9, m.σ, m.dg, m.ϕ, fH2O, (getfield(mechs, :diff)..., x_ϕ_c = getfield(x_ϕ_c_vec, :diff)))
    η_diff = m.σ * 1f9 / ϵ_rate_diff

    return η_diff
end

function get_η_diff(m, viscous_type::Val{HK2003}, params_viscous)
    @unpack mechs, p_dep_calc, ch2o_o, melt_enhancement = params_viscous

    P = p_dep_calc * m.P
    x_ϕ_c_vec = get_melt_settings_for_x_ϕ_c(Val{melt_enhancement}())
    fH2O = calc_fH2O(m.Ch2o, ch2o_o, P, m.T)
    ϵ_rate_diff = sr_flow_law_calculation(m.T, P* 1f9, m.σ, m.dg, m.ϕ, fH2O, (getfield(mechs, :diff)..., x_ϕ_c = getfield(x_ϕ_c_vec, :diff)))
    η_diff = m.σ * 1f9 / ϵ_rate_diff

    return η_diff
end

function get_η_diff(m, viscous_type::Val{xfit_premelt}, params_viscous)
    
    resp_xfit_premelt = forward(
        xfit_premelt(m.T, m.P, m.dg, m.σ, m.ϕ, m.Ch2o, m.T_solidus), params = params_viscous
    )
    
    return resp_xfit_premelt.η
    # requires T_solidus :)
end

function calc_maxwell_times(Gu, m::eburgers_psp, params_btype, JF10_visc, params_viscous, viscous_type, melt_enhancement)

    @unpack TR, PR, dR, E, Vstar, Tau_LR, Tau_HR, Tau_MR, Tau_PR, m_a, m_v, melt_alpha, ϕ_c = params_btype
    x_ϕ_c = getfield(get_melt_settings_for_x_ϕ_c(Val{melt_enhancement}()), :diff)

    if JF10_visc
        scale = (m.dg/ dR)^m_v * exp(E/(MT.gas_R * 1f3) * (1/m.T - 1/TR) + Vstar/(MT.gas_R * 1f3) * (m.P/m.T - PR/TR) * 1f9)
        new_scale = add_melt_effects(m.ϕ, scale, melt_alpha, ϕ_c, x_ϕ_c)
        τ_maxwell = Tau_MR * add_melt_effects(m.ϕ, scale, melt_alpha, ϕ_c, x_ϕ_c)

    else
    # requires η_diff here
        η_diff = get_η_diff(m, Val{viscous_type}(), params_viscous)
        τ_maxwell = η_diff/Gu
    end

    LHP = (m.dg/ dR)^m_a * exp(E/(MT.gas_R * 1f3) * (1/m.T - 1/TR) + Vstar/(MT.gas_R * 1f3) * (m.P/m.T - PR/TR) * 1f9)
    scale_LHP = add_melt_effects(m.ϕ, LHP, melt_alpha, ϕ_c, x_ϕ_c)

    τ_L = Tau_LR * scale_LHP
    τ_H = Tau_HR * scale_LHP
    τ_P = Tau_PR * scale_LHP

    return τ_maxwell, τ_L, τ_H, τ_P
end

function integrate_fn(fn, a,b,N,::Val{:midpoint}) # Defining function for integrating using mid-point rule 
    dx= (b-a)/N;
    mid_points= range(start= a+dx/2, stop= b-dx/2, step= dx); # Mid-points of the intervals 
    f_vals= fn.(mid_points); # function value at the mid-points
    I= dx* sum(f_vals);
    return I;
end

function integrate_fn(fn, a,b,N, ::Val{:trapezoidal}) # Defining function for trapezoid rule 
    dx= (b-a)/N;
    points= range(start= a, stop= b, step= dx); # Edges of the interval
    f_vals= fn.(points); # function value at the edges
    I= 0;
    for i in 1:N # Since the number of points will be N+1, going from 1 to N wil include the first and last points only once
        I+= (f_vals[i]+ f_vals[i+1])/2 
    end
        I*= dx;
    return I; 
end

function integrate_fn(fn, a,b,N, ::Val{:simpson}) # Defining the function for Simpson's rule dx= (b-a)/N;
    points= range(start= a, stop= b, step= dx/2); # Getting the midpoints and the edges of all the intervals 
    f_vals= fn.(points);
    I= 0;
    for i in 0:N-1 # loop over all the intervals, the first interval is defined by i=0, second by i=1, and so on.
        j= 2i+1 # Get the index position of the leftmost point of every interval in the array of all the points (defined as 'points')
        I+= (f_vals[j]+ 4f_vals[j+1]+ f_vals[j+2])/6 
    end
        I*= dx;
    return I; 
end

function integrate_fn(fn, a, b, N, ::Val{:quadgk}) # Defining the function for Simpson's rule dx= (b-a)/N;
    return first(quadgk(fn, a, b))
end

function integrate(J_int_fn, ω::T, integration_params) where {T}
    @unpack l, h, N, rule = integration_params
    # @show rule, typeof
    f(x) = J_int_fn(x, ω)
    return integrate_fn(f, l, h, N, Val{rule}())
end 

function integrate(J_int_fn, ω::T, integration_params) where {T <: AbstractVector{<: Any}}
    return [integrate(J_int_fn, iω, integration_params) for iω in ω]
end

function forward(m::eburgers_psp; params = params_eburgers_psp)
    @unpack integration_params, elastic_type, params_elastic, params_btype, viscous_type, params_viscous, JF10_visc, melt_enhancement = params
    @unpack alf, DeltaB, DeltaP, sig = params_btype

    ω = 2π .* m.f

    resp_elastic = forward(elastic_type([getfield(m, k) for k in fieldnames(elastic_type)]...), params = params_elastic)
    # resp_elastic = forward(anharmonic(m.T, m.P, m.ρ))
    @unpack G, K, Vp, Vs = resp_elastic

    @show G

    Ju = inv(G)

    τ_maxwell, τ_L, τ_H, τ_P = calc_maxwell_times(G, m, params_btype, JF10_visc, params_viscous, Symbol(viscous_type), melt_enhancement)

    J1_int_fn(x, ω) = x^(alf - 1)/(1+(ω*x)^2)

    J2_int_fn(x, ω) = x^alf/(1+(ω*x)^2)

    int_params = (l = τ_L, h = τ_H, N = integration_params.τ_integration_points, rule = integration_params.integration_method) 

    τ_fac = alf * DeltaB/(τ_H^alf - τ_L^alf)

    int1 = τ_fac .* integrate(J1_int_fn, ω, int_params)
    int2 = τ_fac .*ω .* integrate(J2_int_fn, ω, int_params)

    if typeof(ω) <: Vector{<:Any}
        int1 = reshape(int1, size(ω)...)
        int2 = reshape(int2, size(ω)...)
    end

    J1 = 1 .+ int1
    J2 = int2 .+ inv.(ω .* τ_maxwell)

    if DeltaP > 0
        J1_int_fn2(x, ω) = inv(x) * (exp(- log(x/τ_P) * inv(sig) )^2) * 0.5f0 * inv(1 + (ω*x)^2)
        int1 = integrate(J1_int_fn2, ω, (l = 0f0, h= Inf, N = 1, rule = Val{:quadgk}()))
        J1 .= J1 .+ DeltaP .* int1 .* inv.(sig * sqrt(2π))

        J2_int_fn2(x, ω) = (exp(- log(x/τ_P) * inv(sig) )^2) * 0.5f0 * inv(1 + (ω*x)^2)
        int2 = integrate(J2_int_fn2, ω, (l = 0f0, h= Inf, N = 1, rule = Val{:quadgk}()))
        J2 .= J2 .+ DeltaP .* ω * int2 .* inv.(sig * sqrt(2π))
    end

    if typeof(ω) <: Vector{<:Any}
        rmul!(J1, Ju)
        rmul!(J2, Ju)
    end

    J1 = Ju * J1
    J2 = Ju * J2

    Qinv = J2 .* inv.(J1)
    Ma = sqrt.(inv.(J1.^2 + J2.^2))
    Va = sqrt.(Ma ./ m.ρ)

    Vave = sum(Va) * inv(length(m.f))

    return RockPhyAnelastic(
        J1, J2, Qinv, Ma, Va, Vave
    )

end



mutable struct premelt_anelastic{T1, T2, T3, T4, T5, T6, T7, T8, T9}
    T::T1
    P::T2
    dg::T3
    σ::T4
    ϕ::T5
    ρ::T6
    Ch2o::T7
    T_solidus::T8
    f::T9
end

premelt_anelastic(T, P, dg, σ, ϕ, ρ, f) = premelt_anelastic(T, P, dg, σ, ϕ, ρ, 0f0, 0f0, f) # change into kwargs?


function calc_ApSigp(Tn, ϕ, params)
    @unpack α_B, A_B, τ_pp, A_p_fac_1, A_p_fac_2, A_p_fac_3, σ_p_fac_1, σ_p_fac_2, σ_p_fac_3, A_p_Tn_pts, σ_p_Tn_pts, include_direct_melt_effect, β, β_B, poro_Λ = params

    β_p = (include_direct_melt_effect) ? β : 0f0
    
    A_p = 0f0
    if Tn >= A_p_Tn_pts[3]
        A_p = A_p_fac_3 + β_p * ϕ
    elseif Tn < A_p_Tn_pts[3]
        A_p = A_p_fac_3
    elseif Tn < A_p_Tn_pts[2]
        A_p = A_p_fac_1 + A_p_fac_2 * (Tn - A_p_Tn_pts[1])
    elseif Tn < A_p_Tn_pts[1]
        A_p = A_p_fac_1
    end

    σ_p = 0f0
    
    if σ_p_fac_1 < σ_p_Tn_pts[1]
        σ_p = σ_p_fac_1
    elseif (σ_p_fac_1 >= σ_p_Tn_pts[1]) && (σ_p_fac_1 < σ_p_Tn_pts[2])
        σ_p = σ_p_fac_1 + σ_p_fac_2 * (Tn - σ_p_Tn_pts[1])
    else
        σ_p = σ_p_fac_3
    end

    return A_p, σ_p
    
end

function forward(m::premelt_anelastic; params = params_premelt_anelastic)
    # @unpack β2, β2_fit2, α2, τ_cutoff, τ_cutoff_fit2, β1, α_a, α_b, α_c, α_τn, melt_alpha, ϕ_c, x_ϕ_c, scaling_method1,
    
    # @unpack α_B, A_B, τ_pp, A_p_fac_1, A_p_fac_2, A_p_fac_3, σ_p_fac_1, σ_p_fac_2, σ_p_fac_3, A_p_Tn_pts, σ_p_Tn_pts, include_direct_melt_effect, β, β_B, poro_Λ,
    @unpack params_xfit, elastic_type, elastic_params, viscous_params = params

    @unpack include_direct_melt_effect, β_B, poro_Λ, α_B, A_B, τ_pp = params_xfit

    resp_elastic = forward(elastic_type([getfield(m, k) for k in fieldnames(elastic_type)]...), params = elastic_params)
    # resp_elastic = forward(anharmonic(m.T, m.P, m.ρ))
    @unpack G, K, Vp, Vs = resp_elastic

    Ju = inv(G)

    Tn = m.T/ m.T_solidus

    # we only use xfit_premelt vicosity?

    viscous_type = xfit_premelt
    resp_viscous = forward(viscous_type([getfield(m, k) for k in fieldnames(viscous_type)]...), params = viscous_params)
    @unpack ϵ_rate, η = resp_viscous

    τ_m = η * Ju

    A_p, σ_p = calc_ApSigp(Tn, m.ϕ, params_xfit)

    β_B = (include_direct_melt_effect) ? β_B * m.ϕ : 0f0
    
    poro_elastic_factor = (include_direct_melt_effect) ? poro_Λ * m.ϕ : 0f0

    k_temp = A_B + β_B

    pifac = sqrt(π/2)

    T = inv.(m.f)
    p_p = T ./(2π *τ_m)
    ABppa = k_temp .* (p_p).^(α_B)

    lntauapp = log(τ_pp ./ p_p)

    J1 = Ju * (1 + poro_elastic_factor .+ ABppa ./α_B .+ pifac * A_p * σ_p * erfc.(lntauapp ./(sqrt(2f0) * σ_p)))

    J2 = Ju * π/2 * (ABppa .+ A_p * exp(-((lntauapp.^2) ./(2 * σ_p^2))))

    Qinv = J2 .* inv.(J1)
    Ma = sqrt.(inv.(J1.^2 + J2.^2))
    Va = sqrt.(Ma ./ m.ρ)

    Vave = sum(Va) * inv(length(m.f))

    return RockPhyAnelastic(
        J1, J2, Qinv, Ma, Va, Vave
    )
end


m_andrade_psp = andrade_psp(1273f0, 0.2f0, 3.1f0, 10f-3, 0.01f0, 3300f0, 1f0)

forward(m_andrade_psp)

m_eburgers_psp = eburgers_psp(1273f0, 0.2f0, 3.1f0, 10f-3, 0.01f0, 3300f0, 0, 1000f0, 1f0)
# T, P, dg, σ, ϕ, ρ, 0f0, 0f0, f
forward(m_eburgers_psp)


m_premelt_anelastic = premelt_anelastic(1273f0, 0.2f0, 3.1f0, 10f-3, 0.01f0, 3300f0, 0, 1000f0, 1f0)
# T, P, dg, σ, ϕ, ρ, 0f0, 0f0, f
forward(m_premelt_anelastic)

