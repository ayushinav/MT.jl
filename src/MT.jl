module MT
using Optim, LinearAlgebra
using Plots

mutable struct model{T <: Union{AbstractVector{Float32}, AbstractVector{Float64}}}
    m::T
    h::T
end

mutable struct response{T <: Union{AbstractVector{Float32}, AbstractVector{Float64}}}
    ρₐ::T
    ϕ::T
end

# function response(ω::Vector{T}) where T <: Union{Float32, Float64}
#     return response(
#         zero(ω),
#         zero(ω)
#     )
# end

include("fwd_1d.jl")
include("plots.jl")
include("inv.jl")

# global ut::utils
export μ, PRECISION
export model, response
export plot_response, prepare_plot, prepare_plot!
export plot_model, plot_model!
export get_Z, get_appres, get_phase, forward!, forward
export jacobian!, gradient!, lls
# Write your package code here.

end
