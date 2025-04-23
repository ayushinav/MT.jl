"""
create a model type for a given resistivity distribution that can be used to calculate forward response for 1d MT
"""
mutable struct MTModel{T1 <: AbstractArray{<:Any}, T2 <: AbstractArray{<:Any}} <:
               AbstractGeophyModel
    m::T1
    h::T2
end

# explain in blog the reasoning behind this! This covers all 1D, 2D, 3D models for MT.
# make pretty tables to print these models

default_params(::Val{T}) where T <: MT.AbstractGeophyModel = (;)