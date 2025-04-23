struct RockphyElasticDistribution{
    T1 <: Union{Function, Nothing}, T2 <: Union{Function, Nothing}, 
    T3 <: Union{Function, Nothing}, T4 <: Union{Function, Nothing}} <: AbstractRockphyResponseDistribution
    
    G::T1
    K::T2
    Vp::T3
    Vs::T4
end

mutable struct anharmonicDistribution{
    T1 <: Union{Distribution, AbstractArray}, T2 <: Union{Distribution, AbstractArray},
    T3 <: Union{Distribution, AbstractArray}} <: AbstractRockphyModelDistribution
    
    T::T1
    P::T2
    ρ::T3
end

mutable struct anharmonic_poroDistribution{
    T1 <: Union{Distribution, AbstractArray}, T2 <: Union{Distribution, AbstractArray},
    T3 <: Union{Distribution, AbstractArray}, T4 <: Union{Distribution, AbstractArray}} <: AbstractRockphyModelDistribution
    
    T::T1
    P::T2
    ρ::T3
    ϕ::T4
end

mutable struct SLB2005Distribution{
    T1 <: Union{Distribution, AbstractArray}, T2 <: Union{Distribution, AbstractArray}} <: AbstractRockphyModelDistribution
    
    T::T1
    P::T2
end