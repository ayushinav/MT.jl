"""
    mutable struct modelDistribution{T1<: Union{Distribution, AbstractArray}, T2<: Union{Distribution, AbstractArray}} # where T1,T2 
        m::T1
        h::T2
    end

create a placeholder to store the `Distributions.jl` sampler for a priori
"""
mutable struct MTModelDistribution{T1 <: Union{Distribution, AbstractArray},
                                   T2 <: Union{Distribution, AbstractArray}} <:
               AbstractGeophyModelDistribution
    m::T1
    h::T2
end

# we don't need these distributions as structs but as named tuples?
"""
    struct responseDistribution{T1<: Union{Function, Nothing}, T2<: Union{Function, Nothing}} # where T1,T2
        ρₐ::T1
        ϕ::T2
    end

create a placeholder to store functions to obtain `Distributions.jl` samplers for the likelihood function
"""
struct MTResponseDistribution{T1 <: Union{Function, Nothing}, T2 <:
                                                              Union{Function, Nothing}} <:
       AbstractGeophyResponseDistribution
    ρₐ::T1
    ϕ::T2
end

mutable struct RockphyModelDistribution{T1 <: Union{Distribution, AbstractArray},
    T2 <: Union{Distribution, AbstractArray}} <: AbstractRockphyModelDistribution
    params::T1 # vector of parameters 
    p_names::Vector{<:Symbol} # Vector of symbols telling the parameters in vector 
    ϕ::T2 # phase ratios
    model_list::Vector{<:Type}
    mixing_type::Type
    # add water and partition ratios

    # @assert ϕ and model_lists have same length
    # @assert p contains all the variables required by all models in model_list
end

struct RockphyResponseDistribution{T} <: AbstractResponseDistribution
    σ::T
end