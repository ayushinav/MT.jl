abstract type AbstractElasticModel <: AbstractRockphyModel end

## response

mutable struct RockphyElastic{T1,T2,T3,T4} <: AbstractRockphyResponse
    G::T1
    K::T2
    Vp::T3
    Vs::T4
end

"""
T : Temperature (Kelvin)
P : Pressure (GPa)
ρ : density (kg/m³)
"""
mutable struct anharmonic{T1,T2,T3} <: AbstractElasticModel
    T::T1
    P::T2
    ρ::T3
end

mutable struct anharmonic_poro{T1,T2,T3,T4} <: AbstractElasticModel
    T::T1
    P::T2
    ρ::T3
    ϕ::T4
end

mutable struct SLB2005{T1,T2} <: AbstractElasticModel
    T::T1
    P::T2
end

# m_anharmonic = anharmonic(1273f0, 0.2f0, 3300f0)

# forward(m_anharmonic)

# m_anharmonic_poro = anharmonic_poro(1273f0, 0.2f0, 3300f0, 0.01f0)

# forward(m_anharmonic_poro)

# m_SLB2005 = SLB2005_(1273f0, 0.2f0)

# forward(m_SLB2005)

