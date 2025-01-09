"""
`∂(n)`: returns a `n`x`n` matrix for a 1D finite difference stencil using 2 points.
"""
function ∂(n)
    D = I(n) .+ 0.01
    D .= D .- 0.01
    for i in 2:n
        D[i, i - 1] = -1
    end
    return D
end
"""
`χ²(dcal::T, dobs::T; W)`: returns a chi-squared error between the observed and the calculated data. `W` can optionally be passed to weigh points differently.
"""
function χ²(dcal::T, dobs::T;
        W::AbstractMatrix) where {T <: Union{AbstractVector{<:Any}, AbstractVector{<:Any}}}
    sqrt((dcal .- dobs)' * W * (dcal .- dobs) / length(dcal))
end

"""
`struct linear_utils`:
contains the utilities for linearizing the forward model => `mₖ`:model, `Fₖ`: Forward response at `mₖ`, `Jₖ`: Jacobian at `mₖ`.
"""
mutable struct linear_utils{
    T1, T2 <: Union{AbstractVector{Float32}, AbstractVector{Float64}},
    T3 <: Union{AbstractMatrix{Float32}, AbstractMatrix{Float64}}}
    mₖ::T1
    Fₖ::T2
    Jₖ::T3
end

"""
`struct inverse_utils`:
contains the utilities for inversion, once initialized, will not be updated in the inversion iterations `D`: second derivative operator, `W`: weight matrix, `dobs`: data response to be inverted for.
"""
mutable struct inverse_utils{
    T1, T2 <: Union{AbstractMatrix{Float32}, AbstractMatrix{Float64}},
    T3 <: Union{AbstractVector{Float32}, AbstractVector{Float64}}}
    D::T1
    W::T2
    dobs::T3
end

"""
`struct return_code`:
contains the information if the inversion was successful
"""
mutable struct return_code{T1 <: AbstractModel}
    if_pass::Bool
    parameters::NamedTuple
    model_estimate::T1
    misfit_threshold::AbstractFloat
    misfit_achieved::AbstractFloat
end

# linear_utils and inverse_utils are used only in occam, so we do not touch them for now
