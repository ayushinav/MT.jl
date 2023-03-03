const μ= 4π*1e-7;


mutable struct utils
    kj::AbstractArray{Complex{Float64}, 1}
    Z_j::AbstractArray{Complex{Float64}, 1}
    Z_oj::AbstractArray{Complex{Float64}, 1}
end

macro mt_init(nω)
    quote
        nω= $(esc(nω))
        mt_iinit(nω)
    end
end

function mt_iinit(nω)
        global ut= utils(fill(0. *im, nω), fill(0. *im, nω), fill(0. *im, nω));
end

function update_Z!(Z,ω::AbstractArray{Float64, 1},h::AbstractArray{Float64, 1},ρ::AbstractArray{Float64, 1}, Z_oj::AbstractArray{T,1}, kj::AbstractArray{T,1}, Z_j::AbstractArray{T,1}) where T<: ComplexF64
    fill!(Z, 0. *im);
    broadcast!(*, Z, ω, im, μ);
    broadcast!(*, Z, view(ρ, length(ρ)))
    broadcast!(sqrt, Z, Z);
    for j= length(ρ)-1:-1:1
        broadcast!(*, Z_oj, ω, im, μ);#, view(ρ, j));
        broadcast!(*, Z_oj, view(ρ, j));
        broadcast!(sqrt, Z_oj, Z_oj);
        broadcast!(/, kj, Z_oj, view(ρ, j));
        broadcast!(*, kj, -2, view(h,j), kj);
        @. Z_j= Z_oj*((Z_oj+Z)- (Z_oj-Z)*exp(kj))/((Z_oj+Z)+ (Z_oj-Z)*exp(kj));
        @. Z= Z_j;
    end
end

function update_Z!(Z,ω::AbstractArray{Float64, 1},h::AbstractArray{Float64, 1},ρ::AbstractArray{Float64, 1}, ut::utils)
    fill!(Z, 0. *im);
    broadcast!(*, Z, ω, im, μ);
    broadcast!(*, Z, view(ρ, length(ρ)))
    broadcast!(sqrt, Z, Z);
    for j= length(ρ)-1:-1:1
        broadcast!(*, ut.Z_oj, ω, im, μ);#, view(ρ, j));
        broadcast!(*, ut.Z_oj, view(ρ, j));
        broadcast!(sqrt, ut.Z_oj, ut.Z_oj);
        broadcast!(/, ut.kj, ut.Z_oj, view(ρ, j));
        broadcast!(*, ut.kj, -2, view(h,j), ut.kj);
        @. ut.Z_j= ut.Z_oj*((ut.Z_oj+Z)- (ut.Z_oj-Z)*exp(ut.kj))/((ut.Z_oj+Z)+ (ut.Z_oj-Z)*exp(ut.kj));
        @. Z= ut.Z_j;
    end
end

update_Z!(Z,ω::AbstractArray{Float64, 1},h::AbstractArray{Float64, 1},ρ::AbstractArray{Float64, 1})= 
update_Z!(Z,ω,h,ρ, ut.Z_oj, ut.kj, ut.Z_j);
