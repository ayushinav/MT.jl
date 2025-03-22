const model_names_definition = (
    T="Temperature (K)", ρ="Density (kg/m³)", P="Pressure (GPa)",
    ϕ="Porosity", Ch2o_m="Water concentration in melt (ppm)",
    Ch2o_ol="Water concentration in olivine (ppm)",
    Cco2_m="CO₂ concentration in melt (ppm)")

function Base.show(io::IO, m::model) where {model <: AbstractRockphyModel}
    println("Model : ", typeof(m).name.name)
    for k in propertynames(m)
        println(model_names_definition[k], " : ", getfield(m, k))
    end
end

const resp_names_definition = (
    T="Temperature (K)", σ="log₁₀ conductivity (S/m)", G="Elastic shear modulus (Pa)",
    K="Elastic bulk modulus (Pa)", Vp="Elastic P-wave velocity (m/s)",
    Vs="Elastic S-wave velocity (m/s)", ϵ_rate="Strain rate",
    η="Viscosity (Pa s)", J1="Real part of dynamic compliance (Pa⁻¹)",
    J2="Imaginary part of dynamic compliance (Pa⁻¹)", Qinv="Attenuation", M="Modulus (Pa)",
    V="Anelastic S-wave velocity : (m/s)", Vave="Frequency averaged S-wave velocity (m/s)")

function Base.show(io::IO, m::model) where {model <: AbstractRockphyResponse}
    println("Rock physics Response : ", typeof(m).name.name)
    for k in propertynames(m)
        println(resp_names_definition[k], " : ", getfield(m, k))
    end
end

function Base.show(io::IO, m::model) where {model <: model_multiphase}
    if typeof(first(m.mixing_type)) <: MT.single_phase
        println("# Single phase composition")
        println(first(m.model_list))
    elseif typeof(first(m.mixing_type)) <: Union{MT.HS1962_plus, MT.HS1962_minus}
        println("# Two phase composition using ", first(m.mixing_type))
        # for i in 1:2
        println(1 - first(m.ϕ), " : ", m.model_list[1])
        println(first(m.ϕ), " : ", m.model_list[2])
        # end
    end
    println("\n #Parameters")
    for i in eachindex(m.p_names)
        println(m.p_names[i], " : ", m.params[i])
    end
end
