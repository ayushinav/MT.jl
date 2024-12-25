const p_names_definition = (
    T = "Temperature (K)",
    Ch2o_m = "Water concentration in melt (ppm)",
    Ch2o_ol = "Water concentration in olivine (ppm)",
    Cco2_m = "CO₂ concentration in melt (ppm)"
)

function Base.show(io::IO, m::model) where {model <: Union{MT.AbstractMineralModel, MT.AbstractMeltModel}}
    println("Model : ", typeof(m).name.name)
    for k in propertynames(m)
        println(p_names_definition[k], " : ", getfield(m, k))
    end
end

function Base.show(io::IO, m::model) where {model <: mixing_models}
    if typeof(first(m.mixing_type)) <: MT.single_phase
        println("# Single phase composition")
        println(first(m.model_list))
    elseif typeof(first(m.mixing_type)) <: Union{MT.HS1962_plus, MT.HS1962_minus}
        println("# Two phase composition using ", m.mixing_type)
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