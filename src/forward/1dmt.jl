const global μ = 4π * 1.0f-7; # Float32 will promote to Float64 without a problem

"""
`get_Z(ρ,h,ω)`:

returns a tuple of ρₐ and ϕ, given arrays of resistivity `ρ` and thickness `h` for the angular frequenciy `ω`.
"""

pow10(x::T) where {T} = 10^x
const k4_ = (ρₐ=lin_tf2.tf, ϕ=lin_tf2.tf)

function get_Z(ρ::T1, h::T2, ω::T) where {T1, T2, T}
    Z = complex(zero(eltype(ρ)))
    broadcast!(pow10, ρ, ρ)
    # ρ .= 10 .^ ρ
    k = sqrt(im * ω * μ / ρ[end])
    Z = ω * μ / k

    j = length(h)
    @inbounds while j >= 1
        k = sqrt(im * ω * μ / ρ[j])
        Z = ω * μ / k * coth(-im * k * h[j] + acoth(Z / (ω * μ / k)))
        j -= 1
    end

    broadcast!(log10, ρ, ρ)
    Z = conj(Z)
    return get_appres(Z, ω), get_phase(Z)
end

# dispatch on forward for 1d model
"""
`forward(m::model, ω::Vector{T}) where T <: Union{Float32, Float64}`:

returns a  `response` for the given model `m` at the frequencies  `ω`
"""
function forward(m::Tm, ω::T3, response_trans_utils::T) where {Tm <: MTModel, T, T3}
    if !(length(m.h) == length(m.m) - 1)
        error("number of model layers should be 1 less than the number of model parameters")
    end
    n = length(ω)
    ρₐ = zeros(eltype(m.m), n)
    ϕ = zeros(eltype(m.m), n)
    i = 1
    @inbounds while i <= n
        ρₐ[i], ϕ[i] = get_Z2(m.m, m.h, ω[i])
        i += 1
    end
    f1 = response_trans_utils.ρₐ
    f2 = response_trans_utils.ϕ
    MTResponse{typeof(ρₐ), typeof(ϕ)}((f1.(ρₐ)), f2.(ϕ))
end

function forward(m::Tm, ω::T3) where {Tm <: MTModel, T3}
    return forward(m, ω, k4_)
end

# dispatch on forward! for 1d model

"""
`forward!(r::response, m::model, ω::Vector{T}) where T <: Union{Float32, Float64}`:

updates response `r` type for the given model `m` at the frequencies  `ω`
"""
function forward!(r::Tr, m::Tm, ω::T3,
        response_trans_utils::T) where {Tr <: MTResponse, Tm <: MTModel, T, T3}
    if !(length(m.h) == length(m.m) - 1)
        error("number of model layers should be 1 less than the number of model parameters")
    end
    n = length(ω)
    i = 1
    @inbounds while i <= n
        r.ρₐ[i], r.ϕ[i] = get_Z(m.m, m.h, ω[i])
        i += 1
    end
    f1 = response_trans_utils[:ρₐ]
    f2 = response_trans_utils[:ϕ]
    broadcast!(f1, r.ρₐ, r.ρₐ)
    broadcast!(f2, r.ϕ, r.ϕ)
    nothing
end

function forward!(r::Tr, m::Tm, ω::T) where {Tr <: MTResponse, Tm <: MTModel, T}
    forward!(r, m, ω, k4_)
    nothing
end

# the following are defined on scalars, there in-place variants don't make sense

"""
`get_phase(Z)`: returns the phase for impedance
"""
get_phase(Z) = 180 / π * atan(imag(Z) / real(Z));

"""
`get_appres(Z, ω)`: returns the ρₐ for impedance
"""
get_appres(Z, ω) = abs(Z)^2 / (eltype(ω))(μ) / ω;
