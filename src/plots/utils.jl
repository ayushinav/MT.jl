function get_units(m::MTResponse)
    return (ρₐ=:Ωm, ϕ=:ᵒ)
end

function get_scale(m::MTResponse)
    return (ρₐ=:log10, ϕ=:identity)
end

function get_model_labels(m::MTModel) # can be defined for different types of models of other surveys
    return "ρ (Ωm)", "Depth (m)"
end