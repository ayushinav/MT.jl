mutable struct modelDistribution{T1<: Union{Distribution, AbstractArray}, T2<: Union{Distribution, AbstractArray}} # where T1,T2 
    m::T1
    h::T2
end

# we don't need these distributions as structs but as named tuples?
struct responseDistribution{T1<: Union{Function, Nothing}, T2<: Union{Function, Nothing}} # where T1,T2
    ρₐ::T1
    ϕ::T2
end