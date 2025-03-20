
mutable struct RockphyViscous{T1,T2} #<: AbstractRockphyResponse
    ϵ_rate::T1
    η::T2
end

"""
T : K
P : GPa
dg : μm
σ : GPa
ϕ
"""
mutable struct HZK2011{T1, T2, T3, T4, T5}
    T::T1
    P::T2
    dg::T3
    σ::T4
    ϕ::T5
end


"""
T : K
P : GPa
dg : μm
σ : GPa
ϕ
"""
mutable struct HK2003{T1, T2, T3, T4, T5, T6}
    T::T1
    P::T2
    dg::T3
    σ::T4
    ϕ::T5
    Ch2o::T6
end

HK2003(T, P, dg, σ, ϕ) = HK2003(T, P, dg, σ, ϕ, 0f0)

# TODO : figure this Ch2o out


mutable struct xfit_premelt{T1, T2, T3, T4, T5, T6, T7}
    T::T1
    P::T2
    dg::T3
    σ::T4
    ϕ::T5
    Ch2o::T6
    T_solidus::T7
end

xfit_premelt(T, P, dg, σ, ϕ, T_solidus) = xfit_premelt(T, P, dg, σ, ϕ, 0f0, T_solidus)




# m_HZK2011 = HZK2011(1273f0, 0.2f0, 3.1f0, 10f-3, 0.01f0)

# forward(m_HZK2011)

# m_HK2003 = HK2003(1273f0, 0.2f0, 3.1f0, 10f-3, 0.01f0)

# forward(m_HK2003)

# m_xfit_premelt = xfit_premelt(1273f0, 0.2f0, 3.1f0, 10f-3, 0.01f0, 1000f0)

# forward(m_xfit_premelt)

# # RockphyViscous{Float32, Float64}(0.0f0, 7.133103851381647e13)