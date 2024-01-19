module MT
using Optim, LinearAlgebra
import Plots: plot, plot!
mutable struct model
    ρ::Vector{Float64}
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
export plot, plot!
export get_Z, get_appres, get_phase, forward!, forward
export jacobian!, gradient!, lls
# Write your package code here.

end
