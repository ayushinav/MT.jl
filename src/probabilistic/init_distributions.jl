mutable struct modelDistribution{T1<: Union{Distribution, AbstractArray}, T2<: Union{Distribution, AbstractArray}} # where T1,T2 
    m::T1
    h::T2
end

# we don't need these distributions as structs but as named tuples?
mutable struct responseDistribution{T1<: Union{Distribution, Nothing}, T2<: Union{Distribution, Nothing}} # where T1,T2
    ρₐ::T1
    ϕ::T2
end