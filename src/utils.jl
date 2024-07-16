# utils to help bump for Abstract types

import LinearAlgebra: zero

function zero(x::resp) where {resp <: AbstractResponse}
    typeof(x)([zero(getfield(x, k)) for k in fieldnames(resp)]...)
end;
function zero(x::model) where {model <: AbstractModel}
    typeof(x)([zero(getfield(x, k)) for k in fieldnames(model)]...)
end;

# zero(x::MTResponse) = MTResponse([zero(getfield(x, k)) for k in fieldnames(MTResponse)]...)
# zero(x::MTModel) = MTModel([zero(getfield(x, k)) for k in fieldnames(MTResponse)]...)

import Base: copy
function copy(x::resp) where {resp <: AbstractResponse}
    typeof(x)([copy(getfield(x, k)) for k in fieldnames(resp)]...)
end;
function copy(x::model) where {model <: AbstractModel}
    typeof(x)([copy(getfield(x, k)) for k in fieldnames(model)]...)
end;

# copy(x::MTResponse) = MTResponse([copy(getfield(x, k)) for k in fieldnames(MTResponse)]...)
# copy(x::MTModel) = MTModel([copy(getfield(x, k)) for k in fieldnames(MTResponse)]...)

#=
 This works for now, dispatching on each model/response type, but can we make it general for all types? Something like
    copy(x<:AbstractGeophyModel) = typeof(x)([copy(getfield(x, k)) for k in fieldnames(MTResponse)]...)
Should work probably, but we'll look into. Also 
=#

# probably only required in occam codes, maybe a generated function is needed?

function zero_abstract(m::mtresponse) where {
                                             mtresponse <:
                                             MTResponse{<:AbstractVector{<:Any},
                                                        <:AbstractVector{<:Any}}}
    MTResponse{AbstractVector, AbstractVector}(zero(m.ρₐ), zero(m.ϕ))
end

# this probably does it for all 1D, 2D, 3D models and their corresponding responses? Atleast 1d for sure, some brainstorming maybe required for 2D and 3D.
# Again, can we make a generated function here?

function inverse(t::mtresponse; abstract=false) where {mtresponse <: MTResponse}
    if abstract
        return MTModel{AbstractArray{<:Any, length(size(t.ρₐ))},
                       AbstractArray{<:Any, length(size(t.ρₐ))}}
    else
        vec_type = typeof(t.ρₐ)
        return MTModel{vec_type, vec_type}
    end
end

function sample(d::mtdist) where {mtdist <: MTModelDistribution}
    vec_type = typeof(rand(d.m))
    return MTModel{vec_type, vec_type}
end

function forward(t::mtmodel; abstract=false) where {mtmodel <: MTModel} # {<:AbstractVector{<:Any}, <:AbstractVector{<:Any}}}
    if abstract
        return MTResponse{AbstractArray{<:Any, length(sie(m.m))},
                          AbstractArray{<:Any, length(sie(m.m))}}
    else
        vec_type = typeof(t.m)
        return MTResponse{vec_type, vec_type}
    end
end