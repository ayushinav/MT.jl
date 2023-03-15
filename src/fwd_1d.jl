const μ= 4π*1e-7;


mutable struct utils
    k::AbstractArray{Complex{Float64}, 1}
    # Ztmp::AbstractArray{Complex{Float64}, 1}
    # Z_oj::AbstractArray{Complex{Float64}, 1}
end

macro mt_init(nω)
    quote
        nω= $(esc(nω))
        mt_iinit(nω)
    end
end

function mt_iinit(nω)
        global ut= utils(fill(0. *im, nω));
end


## iip function to get impedance
function get_Z!(Z,ω::AbstractArray{Float64, 1},h::AbstractArray{Float64, 1},ρ::AbstractArray{Float64, 1}, ut::utils)
    
    # k= @SVector fill(im.*0, size(ω)...)
    @. ut.k= sqrt(im*ω*μ/ρ[end])
    
    fill!(Z, im.*0);
    @. Z= ω*μ/ut.k;

    @assert length(h)== length(ρ)- 1
    
    for j in length(h):-1:1
        @. ut.k= sqrt(im*ω*μ/ρ[j])
        @. Z= ω*μ/ut.k*coth(-im*ut.k*h[j] + acoth(Z/(ω*μ/ut.k)))
    end
    conj!(Z);
    return nothing;
end

get_Z!(Z,ω::AbstractArray{Float64, 1},h::AbstractArray{Float64, 1},ρ::AbstractArray{Float64, 1})=
get_Z!(Z, ω, h, ρ, ut)
    

## oop function to get impedance
function get_Z(ω, h, ρ)
    Z= fill(0. *im, size(ω));
    get_Z!(Z, ω, ρ, h);
    return Z;
end
