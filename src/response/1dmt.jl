mutable struct MTResponse{T1, T2} <:
               AbstractGeophyResponse
    ρₐ::T1
    ϕ::T2
end

#=
More often than not, T1 and T2 will be same, that is typeof(ρₐ) and typeof(ϕ) will be same but this might become different when using a rock physics model in front of this. 
For now, we go with the same design as for `MTModel`.
=#
