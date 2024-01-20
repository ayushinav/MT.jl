module MT
using Optim, LinearAlgebra
using Plots
mutable struct model
    m::Vector{Float64}
    h::Vector{Float64}
end

mutable struct response
    ρₐ::Vector{Float64}
    ϕ::Vector{Float64}
end

function response(ω::Vector{T}) where T <: Union{Float32, Float64}
    return response(
        zero(ω),
        zero(ω)
    )
end

include("fwd_1d.jl")
include("plots.jl")
include("inv.jl")

# global ut::utils
export μ
export model, response
export plot_response, prepare_plot, prepare_plot!
export plot_model, plot_model!
export get_Z, get_appres, get_phase, forward!, forward
export jacobian!, gradient!, lls
# Write your package code here.

end
