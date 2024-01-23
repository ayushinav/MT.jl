global const μ= 4π*1f-7; # Float32 will promote to Float64 without a problem

"""
`get_Z(ρ,h,ω)`:
returns a tuple of ρₐ and ϕ, given arrays of resistivity `ρ` and thickness `h` for the angular frequenciy `ω`.
"""
function get_Z(ρ::AbstractVector{T}, h::AbstractVector{T}, ω::T) where T <: Union{Float32, Float64}
    
    Z= complex(zero(eltype(ρ)));
    k= sqrt(im*ω*μ/ρ[end])
    Z= ω*μ/k;

    j=length(h);
     @inbounds while j>=1
        k= sqrt(im*ω*μ/ρ[j])
        Z= ω*μ/k*coth(-im*k*h[j] + acoth(Z/(ω*μ/k)))
        j-=1;
    end
    Z= conj(Z);
    return get_appres(Z, ω), get_phase(Z);
end
"""
`forward(m::model, ω::Vector{T}) where T <: Union{Float32, Float64}`:
returns a  `response` for the given model `m` at the frequencies  `ω`
"""
function forward(m::model, ω::Vector{T}) where T <: Union{Float32, Float64} # don't maek ω an abstract vector because we do the forward modeling by accessing each element, or maybe we can
    # the following line check is why we do not use the same fn name here, so that the checks happen just once for all the frequencies.
    if !(length(m.h)== length(m.m)- 1)
        error("number of model layers should be 1 less than the number of model parameters")
    end
    n= length(ω);
    ρₐ= zeros(eltype(m.m), n);
    ϕ= zeros(eltype(m.m), n);
    i=1;
    @inbounds while i<=n
        ρₐ[i], ϕ[i]= get_Z(m.m, m.h, ω[i]);
        i+=1;
    end
    return response(ρₐ, ϕ);
end

"""
`forward!(r::response, m::model, ω::Vector{T}) where T <: Union{Float32, Float64}`:
updates response `r` type for the given model `m` at the frequencies  `ω`
"""
function forward!(r::response, m::model, ω::Vector{T}) where T <: Union{Float32, Float64}
    if !(length(m.h)== length(m.m)- 1)
        error("number of model layers should be 1 less than the number of model parameters")
    end
    n= length(ω);
    i=1;
    @inbounds while i<=n
        r.ρₐ[i], r.ϕ[i]= get_Z(m.m, m.h, ω[i]);
        i+=1;
    end
    nothing;
end

# the following are defined on scalars, there in-place variants don't make sense

"""
`get_phase(Z)`: returns the phase for impedance
"""
get_phase(Z)= 180/π* atan(imag(Z)/real(Z));

"""
`get_appres(Z, ω)`: returns the ρₐ for impedance
"""
get_appres(Z, ω)= abs(Z)^2/(eltype(ω))(μ)/ω;