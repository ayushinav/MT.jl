struct RockphyViscousDistribution{
    T1 <: Union{Function, Nothing}, T2 <: Union{Function, Nothing}} <:
       AbstractRockphyResponseDistribution
    ϵ_rate::T1
    η::T2
end

mutable struct HZK2011Distribution{
    T1 <: Union{Distribution, AbstractArray}, T2 <: Union{Distribution, AbstractArray},
    T3 <: Union{Distribution, AbstractArray}, T4 <: Union{Distribution, AbstractArray},
    T5 <: Union{Distribution, AbstractArray}} <: AbstractRockphyModelDistribution
    T::T1
    P::T2
    dg::T3
    σ::T4
    ϕ::T5
end

mutable struct HK2003Distribution{
    T1 <: Union{Distribution, AbstractArray}, T2 <: Union{Distribution, AbstractArray},
    T3 <: Union{Distribution, AbstractArray}, T4 <: Union{Distribution, AbstractArray},
    T5 <: Union{Distribution, AbstractArray}, T6 <: Union{Distribution, AbstractArray}} <:
               AbstractRockphyModelDistribution
    T::T1
    P::T2
    dg::T3
    σ::T4
    ϕ::T5
    Ch2o_ol::T6
end

mutable struct xfit_premeltDistribution{
    T1 <: Union{Distribution, AbstractArray}, T2 <: Union{Distribution, AbstractArray},
    T3 <: Union{Distribution, AbstractArray}, T4 <: Union{Distribution, AbstractArray},
    T5 <: Union{Distribution, AbstractArray}, T6 <: Union{Distribution, AbstractArray}} <:
               AbstractRockphyModelDistribution
    T::T1
    P::T2
    dg::T3
    σ::T4
    ϕ::T5
    T_solidus::T6
end
