# utils to help bump for Abstract types

import LinearAlgebra:zero

zero(x::resp) where {resp <: AbstractResponse} = typeof(x)([zero(getfield(x, k)) for k in fieldnames(resp)]...);
zero(x::model) where {model <: AbstractModel} = typeof(x)([zero(getfield(x, k)) for k in fieldnames(model)]...);

# zero(x::MTResponse) = MTResponse([zero(getfield(x, k)) for k in fieldnames(MTResponse)]...)
# zero(x::MTModel) = MTModel([zero(getfield(x, k)) for k in fieldnames(MTResponse)]...)

import Base:copy
copy(x::resp) where {resp <: AbstractResponse} = typeof(x)([copy(getfield(x, k)) for k in fieldnames(resp)]...);
copy(x::model) where {model <: AbstractModel} = typeof(x)([copy(getfield(x, k)) for k in fieldnames(model)]...);

# copy(x::MTResponse) = MTResponse([copy(getfield(x, k)) for k in fieldnames(MTResponse)]...)
# copy(x::MTModel) = MTModel([copy(getfield(x, k)) for k in fieldnames(MTResponse)]...)


#=
 This works for now, dispatching on each model/response type, but can we make it general for all types? Something like
    copy(x<:AbstractGeophyModel) = typeof(x)([copy(getfield(x, k)) for k in fieldnames(MTResponse)]...)
Should work probably, but we'll look into. Also 
=#


# probably only required in occam codes, maybe a generated function is needed?

function zero_abstract(m::mtresponse) where {mtresponse <: MTResponse{<:AbstractVector{<:Any}, <:AbstractVector{<:Any}}}
    MTResponse{AbstractVector, AbstractVector}(
        zero(m.ρₐ),
        zero(m.ϕ)
    )
end


