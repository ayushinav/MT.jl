"""
create a model type for a given resistivity distribution that can be used to calculate forward response for 1d MT
"""
mutable struct model{T}
    m::T
    h::T
end