mutable struct MTResponse{T1, T2} <: AbstractGeophyResponse
    ρₐ::T1
    ϕ::T2
end

#=
More often than not, T1 and T2 will be same, that is typeof(ρₐ) and typeof(ϕ) will be same but this might become different when using a rock physics model in front of this. 
For now, we go with the same design as for `MTModel`.
=#

"""
create a model type for a given resistivity distribution 
that can be used to calculate forward response for 1d MT

## Usage
```jldoctest
julia> m = [4., 2., 3.];
julia> h = [1000., 100.];

julia> MTModel(m, h)

1D MTModel : 
Layer    log(ρ)  h
____________________________
1        4.0     1000.0
2        2.0     100.0
3        0.477   ∞
```
"""
mutable struct MTModel{T1 <: AbstractArray{<:Any}, T2 <: AbstractArray{<:Any}} <:
               AbstractGeophyModel
    m::T1
    h::T2
end

# explain in blog the reasoning behind this! This covers all 1D, 2D, 3D models for MT.
# make pretty tables to print these models

default_params(::Type{T}) where {T <: MT.AbstractGeophyModel} = (;)
