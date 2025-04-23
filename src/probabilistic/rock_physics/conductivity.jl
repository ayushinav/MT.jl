struct RockphyCondDistribution{T <: Union{Function, Nothing}} <:
       AbstractRockphyResponseDistribution
    σ::T
end

mutable struct SEO3Distribution{F <: Union{Distribution, AbstractArray}} <:
               AbstractRockphyModelDistribution
    T::F
end

mutable struct UHO2014Distribution{
    F1 <: Union{Distribution, AbstractArray}, F2 <: Union{Distribution, AbstractArray}} <:
               AbstractRockphyModelDistribution
    T::F1
    Ch2o_ol::F2
end

mutable struct Jones2012Distribution{
    F1 <: Union{Distribution, AbstractArray}, F2 <: Union{Distribution, AbstractArray}} <:
               AbstractRockphyModelDistribution
    T::F1
    Ch2o_ol::F2
end

mutable struct Poe2010Distribution{
    F1 <: Union{Distribution, AbstractArray}, F2 <: Union{Distribution, AbstractArray}} <:
               AbstractRockphyModelDistribution
    T::F1
    Ch2o_ol::F2
end

mutable struct Wang2006Distribution{
    F1 <: Union{Distribution, AbstractArray}, F2 <: Union{Distribution, AbstractArray}} <:
               AbstractRockphyModelDistribution
    T::F1
    Ch2o_ol::F2
end

mutable struct Yoshino2009Distribution{
    F1 <: Union{Distribution, AbstractArray}, F2 <: Union{Distribution, AbstractArray}} <:
               AbstractRockphyModelDistribution
    T::F1
    Ch2o_ol::F2
end

mutable struct const_matrixDistribution{F <: Union{Distribution, AbstractArray}} <:
               AbstractRockphyModelDistribution
    σ::F
end

mutable struct Ni2011Distribution{
    F1 <: Union{Distribution, AbstractArray}, F2 <: Union{Distribution, AbstractArray}} <:
               AbstractRockphyModelDistribution
    T::F1
    Ch2o_m::F2
end

mutable struct Sifre2014Distribution{
    F1 <: Union{Distribution, AbstractArray}, F2 <: Union{Distribution, AbstractArray},
    F3 <: Union{Distribution, AbstractArray}} <: AbstractRockphyModelDistribution
    T::F1
    Ch2o_m::F2
    Cco2_m::F2
end

mutable struct Gaillard2008Distribution{F <: Union{Distribution, AbstractArray}} <:
               AbstractRockphyModelDistribution
    T::F
end

