const μ= 4π*1e-7;

function get_Z!(Z,ρ::AbstractArray{Float64, 1},h::AbstractArray{Float64, 1},ωs::AbstractArray{Float64, 1})
    @assert length(h)== length(ρ)- 1
    for i in 1:length(ωs)
        k= sqrt(im*ωs[i]*μ/ρ[end])
        Z[i]= ωs[i]*μ/k;
        for j in length(h):-1:1
            k= sqrt(im*ωs[i]*μ/ρ[j])
            Z[i]= ωs[i]*μ/k*coth(-im*k*h[j] + acoth(Z[i]/(ωs[i]*μ/k)))
        end
    end
    conj!(Z);
    return nothing;
end

function get_Z(ρ, h, ωs)
    Z= fill(im.*0., size(ωs));
    @assert length(h)== length(ρ)- 1
    for i in 1:length(ωs)
        k= sqrt(im*ωs[i]*μ/ρ[end])
        Z[i]= ωs[i]*μ/k;
        for j in length(h):-1:1
            k= sqrt(im*ωs[i]*μ/ρ[j])
            Z[i]= ωs[i]*μ/k*coth(-im*k*h[j] + acoth(Z[i]/(ωs[i]*μ/k)))
        end
        Z[i]= conj(Z[i]);
    end

    return Z
end

function get_appres!(ρₐ, Z, ω)
    # abs.(Z_aray)./ω./μ
    for i in 1:length(Z)
        ρₐ[i]= abs(Z[i])^2/μ/ω[i]
    end
    return nothing;
end

function get_appres(Z, ω)
    ρₐ= zeros(size(Z))
    for i in 1:length(Z)
        ρₐ[i]= abs(Z[i])^2/μ/ω[i]
    end
    return ρₐ
end

function get_phase!(ϕ, Z)
    for i in 1:length(Z)
        ϕ[i]= 180/π* atan(imag(Z[i])/real(Z[i]))
    end
    return nothing;
end

function get_phase(Z)
    ϕ= zeros(size(Z))
    for i in 1:length(Z)
        ϕ[i]= 180/π* atan(imag(Z[i])/real(Z[i]))
    end
    return ϕ
end

function forward!(d::response, m::model, ω)
    get_Z!(d.Z, m.ρ, m.h, ω)
    get_appres!(d.ρₐ, d.Z, ω)
    get_phase!(d.ϕ, d.Z)
    return nothing;
end

function forward(m::model, ω)
    d= response(fill(im.*0, length(ω)),
        zeros(length(ω)),
        zeros(length(ω))
    )
    d.Z= get_Z(m.ρ, m.h, ω)
    d.ρₐ= get_appres(d.Z, ω)
    d.ϕ= get_phase(d.Z)
    
    return d;
end