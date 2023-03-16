module MT
import Plots: plot, plot!
struct model
    ρ::Vector{Float64}
    h::Vector{Float64}
end

struct response
    Z::Vector{ComplexF64}
    ρₐ::Vector{Float64}
    ϕ::Vector{Float64}
end
include("fwd_1d.jl")
include("plots.jl")

# global ut::utils
export μ
export model, response
export plot, plot!
export get_Z!, get_Z, get_appres!, get_appres, get_phase!, get_phase!, forward!, forward
# Write your package code here.

end
