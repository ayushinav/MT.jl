struct RockPhyAnelasticDistribution{
    T1 <: Union{Distribution, AbstractArray}, T2 <: Union{Distribution, AbstractArray},
    T3 <: Union{Distribution, AbstractArray}, T4 <: Union{Distribution, AbstractArray},
    T5 <: Union{Distribution, AbstractArray}, T6 <: Union{Distribution, AbstractArray}} <:
       AbstractRockphyResponseDistribution
    J1::T1
    J2::T2
    Qinv::T3
    M::T4
    V::T5
    Vave::T6
end

mutable struct andrade_pspDistribution{
    T1 <: Union{Distribution, AbstractArray}, T2 <: Union{Distribution, AbstractArray},
    T3 <: Union{Distribution, AbstractArray}, T4 <: Union{Distribution, AbstractArray},
    T5 <: Union{Distribution, AbstractArray}, T6 <: Union{Distribution, AbstractArray},
    T7 <: Union{Distribution, AbstractArray}} <: AbstractRockphyModelDistribution
    T::T1
    P::T2
    dg::T3
    σ::T4
    ϕ::T5
    ρ::T6
    f::T7
end

mutable struct eburgers_pspDistribution{
    T1 <: Union{Distribution, AbstractArray}, T2 <: Union{Distribution, AbstractArray},
    T3 <: Union{Distribution, AbstractArray}, T4 <: Union{Distribution, AbstractArray},
    T5 <: Union{Distribution, AbstractArray}, T6 <: Union{Distribution, AbstractArray},
    T7 <: Union{Distribution, AbstractArray}, T8 <: Union{Distribution, AbstractArray},
    T9 <: Union{Distribution, AbstractArray}} <: AbstractRockphyModelDistribution
    T::T1
    P::T2
    dg::T3
    σ::T4
    ϕ::T5
    ρ::T6
    Ch2o_ol::T7
    T_solidus::T8
    f::T9
end

mutable struct premelt_anelasticDistribution{
    T1 <: Union{Distribution, AbstractArray}, T2 <: Union{Distribution, AbstractArray},
    T3 <: Union{Distribution, AbstractArray}, T4 <: Union{Distribution, AbstractArray},
    T5 <: Union{Distribution, AbstractArray}, T6 <: Union{Distribution, AbstractArray},
    T7 <: Union{Distribution, AbstractArray}, T8 <: Union{Distribution, AbstractArray},
    T9 <: Union{Distribution, AbstractArray}} <: AbstractRockphyModelDistribution
    T::T1
    P::T2
    dg::T3
    σ::T4
    ϕ::T5
    ρ::T6
    Ch2o_ol::T7
    T_solidus::T8
    f::T9
end

mutable struct xfit_mxwDistribution{
    T1 <: Union{Distribution, AbstractArray}, T2 <: Union{Distribution, AbstractArray},
    T3 <: Union{Distribution, AbstractArray}, T4 <: Union{Distribution, AbstractArray},
    T5 <: Union{Distribution, AbstractArray}, T6 <: Union{Distribution, AbstractArray},
    T7 <: Union{Distribution, AbstractArray}, T8 <: Union{Distribution, AbstractArray},
    T9 <: Union{Distribution, AbstractArray}} <: AbstractRockphyModelDistribution
    T::T1
    P::T2
    dg::T3
    σ::T4
    ϕ::T5
    ρ::T6
    Ch2o_ol::T7
    T_solidus::T8
    f::T9
end

mutable struct andrade_analyticalDistribution{
    T1 <: Union{Distribution, AbstractArray}, T2 <: Union{Distribution, AbstractArray},
    T3 <: Union{Distribution, AbstractArray}, T4 <: Union{Distribution, AbstractArray},
    T5 <: Union{Distribution, AbstractArray}, T6 <: Union{Distribution, AbstractArray},
    T7 <: Union{Distribution, AbstractArray}, T8 <: Union{Distribution, AbstractArray},
    T9 <: Union{Distribution, AbstractArray}} <: AbstractRockphyModelDistribution
    T::T1
    P::T2
    dg::T3
    σ::T4
    ϕ::T5
    ρ::T6
    Ch2o_ol::T7
    T_solidus::T8
    f::T9
end
