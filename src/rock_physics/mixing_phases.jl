
"""
mixing_models(params, ϕ, [::SEO3, ::Ni2011], ::HS_plus)
"""
# ::HS_plus type doesn't necessarily have to store anything inside
mutable struct mixing_models{T1, T2} <: AbstractRockphyModel
    params::T1 # vector of parameters 
    p_names::Vector{<:Symbol} # Vector of symbols telling the parameters in vector 
    ϕ::T2 # phase ratios
    model_list::Vector{<:Type}
    mixing_type::Type
    # add water and partition ratios

    # @assert ϕ and model_lists have same length
    # @assert p contains all the variables required by all models in model_list
end

# we can have total water conc and the ratio by which it goes into different phases
# this would just imply adding one more argument while calling `forward` for different rock physics types

function forward(m::model, p) where {model <: mixing_models}

    σs = []; #zeros(length(m.model_list))

    params = (; zip(
        m.p_names,
        m.params
    )...)

    # water partition code

    for i in eachindex(m.model_list)
        var_list = [(fieldnames(m.model_list[i]))...];
        σ = forward(m.model_list[i](params[var_list]...), [])
        push!(σs, σ)
    end

    σ_net =  mix_models(σs, m.ϕ, m.mixing_type())

    return RockphyCond([σ_net])

end


function mix_models(σs, ϕ, ::HS_plus)
    σ_max = maximum(σs);
    # @show size(σs), size(ϕ)

    σ_plus = inv(sum(ϕ .* inv.(σ_max .+ σs))) - σ_max
    return σ_plus
end

function mix_models(σs, ϕ, ::HS_minus)
    σ_min = minimum(σs);

    σ_minus = inv(sum(ϕ .* inv.(σ_max .+ σs))) - σ_min
    return σ_minus
end

function mix_models(σs, ϕ, ::single_phase)
    @assert length(σs) == 1
    return first(σs)
end