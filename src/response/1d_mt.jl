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