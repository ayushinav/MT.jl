"""
create a model type for a given resistivity distribution that can be used to calculate forward response for 1d MT
"""
mutable struct model{T <: Union{AbstractVector{Float32}, AbstractVector{Float64}}}
    m::T
    h::T
end