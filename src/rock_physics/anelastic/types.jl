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

mutable struct xfit_mxw{T1, T2, T3, T4, T5, T6, T7, T8, T9}
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

mutable struct andrade_analytical{T1, T2, T3, T4, T5, T6, T7, T8, T9}
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


# am_andrade_psp = andrade_psp(1273f0, 0.2f0, 3.1f0, 10f-3, 0.01f0, 3300f0, 1f0)

# forward(m_andrade_psp)

# m_eburgers_psp = eburgers_psp(1273f0, 0.2f0, 3.1f0, 10f-3, 0.01f0, 3300f0, 0, 1000f0, 1f0)
# # T, P, dg, σ, ϕ, ρ, 0f0, 0f0, f
# forward(m_eburgers_psp)

# m_premelt_anelastic = premelt_anelastic(1273f0, 0.2f0, 3.1f0, 10f-3, 0.01f0, 3300f0, 0, 1000f0, 1f0)
# # T, P, dg, σ, ϕ, ρ, 0f0, 0f0, f
# forward(m_premelt_anelastic)

# m_xfit_mxw = xfit_mxw(1273f0, 0.2f0, 3.1f0, 10f-3, 0.01f0, 3300f0, 1000f0, 0, 1f0)
# # T, P, dg, σ, ϕ, ρ, 0f0, 0f0, f
# forward(m_xfit_mxw)

# m_andrade_analytical = andrade_analytical(1273f0, 0.2f0, 3.1f0, 10f-3, 0.01f0, 3300f0, 1000f0, 0, 1f0)
# # T, P, dg, σ, ϕ, ρ, 0f0, 0f0, f
# forward(m_andrade_analytical)