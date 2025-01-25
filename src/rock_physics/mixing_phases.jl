"""
    construct_mixing_models(params, p_names, ϕ, model_list, mixing_type)

returns a `mixing_models` type containing all the variables for rock physics modeling

## Arguments

  - `params` : parameter values for rock physics models

  - `p_names` : list of symbols, each associating with each value in `params`
  - `ϕ` : vol fraction of different phases

      + for `mixing_type` = `single_phase`, `ϕ` = `[1]` (total vol occupied by one phase)
      + for `mixing_type` = `HS1962_plus` and `HS1962_minus`, `ϕ` is a vector of one element that can be varied but should be less than 1.
        This is the vol melt fraction.
  - `model_list` : list of model types required to build up the rock physics model and mixing them
  - `mixing_type` : determines how to mix model types. Current mixing laws include

      + `single_phase` : no mixing, when there's only a single phase
      + `HS1962_plus` : mixing two phases to get the Hashim Strikman upper bound
      + `HS1962_minus` : mixing two phases to get the Hashim Strikman lower bound
      + `MAL` : mixing the two phases using the Modified Archie's Law
"""
function construct_mixing_models(
        params::Vector, p_names::Vector, ϕ::Vector, model_list::Vector, mixing_type::Vector)
    var_list = vcat([[fieldnames(ir)...] for ir in model_list]...)
    unique!(var_list)
    f_ = reduce(&, [ir ∈ p_names for ir in var_list])
    @assert f_==true """all the variables required by models are not included in `p_names`. \n
    parameters mentioned : $p_names \n
    all parameters required : $var_list \n
    """

    if typeof(first(mixing_type)) <: Union{HS1962_plus, HS1962_minus, MAL}
        @assert length(model_list)==2 "`$mixing_type` model allows for only 2 models to mix, the first one being the solid and the second melt"
        @assert length(ϕ)==1 "only the fraction of the second component, ie melt, $(model_list[2]) is needed"
        @assert length(params)==length(p_names) "mismatch between the number of parameter names and their values"

        return mixing_models(params, p_names, ϕ, model_list, mixing_type)

    elseif typeof(first(mixing_type)) <: single_phase
        @assert length(model_list)==1 "single phase models require only one model in `model_list`"
        @assert length(ϕ) == 1&&first(ϕ) == 1 "single phase models will have total fraction for the single phase"
        return mixing_models(params, p_names, ϕ, model_list, mixing_type)

        return mixing_models(params, p_names, ϕ, model_list, mixing_type)

    else
        @assert length(model_list)==(length(ϕ) + 1) """each phase and it's vol fraction should be provided, except for the last one, 
        where the vol fraction for the last one is obtained from the remaining variables 1 - ∑ϕ"""
        return mixing_models(params, p_names, ϕ, model_list, mixing_type)
    end
end

"""
    mixing_models

constructs a `mixing_models` type which can then be used to do rock physics modeling. Should be called using `construct_mixing_models`
"""
mutable struct mixing_models{T1, T2} <: AbstractRockphyModel
    params::T1 # vector of parameters 
    p_names::Vector{Symbol} # Vector of symbols telling the parameters in vector
    ϕ::T2 # phase ratios
    model_list::Vector{<:Any}
    mixing_type
    # Ch2o::T4
    # h2o_part::T5
end

# we can have total water conc and the ratio by which it goes into different phases
# this would just imply adding one more argument while calling `forward` for different rock physics types

function forward(m::model, p) where {model <: mixing_models}
    σs = [] #zeros(length(m.model_list))

    params = (; zip(m.p_names, m.params)...)

    # water partition code

    for i in eachindex(m.model_list)
        var_list = [(fieldnames(m.model_list[i]))...]
        σ = forward(m.model_list[i](params[var_list]...), [])
        push!(σs, σ)
    end

    σ_net = mix_models(σs, m.ϕ, first(m.mixing_type))

    return RockphyCond([σ_net])
end

function mix_models(σs, ϕ, ::HS1962_plus)
    σ_max = 10.0f0^maximum(σs)
    σ_min = 10.0f0^minimum(σs)
    phi = first(ϕ)

    num = 3 * (1 - phi) * (σ_max - σ_min) # numerator
    den = 3 * σ_max - phi * (σ_max - σ_min) # denominator
    esig = σ_max * (1 - (num / den))

    return log10(esig)
end

function mix_models(σs, ϕ, ::HS1962_minus)
    σ_max = 10.0f0^maximum(σs)
    σ_min = 10.0f0^minimum(σs)
    phi = first(ϕ)

    num = 3 * (phi) * (σ_max - σ_min) # numerator
    den = 3 * σ_min + (1 - phi) * (σ_max - σ_min) # denominator
    esig = σ_min * (1 + (num / den))

    return log10(esig)
end

function mix_models(σs, ϕ, mal::MAL)
    σ_fluid = 10.0f0^(σs[2])
    σ_matrix = 10.0f0^(σs[1])

    phi = first(ϕ)
    sig = σ_fluid

    if phi < 1
        p = log10(1 - phi^mal.m) * inv(log10(1 - phi))
        sig = σ_fluid * phi^mal.m + σ_matrix * (1 - phi)^p
    end

    return log10(sig)
end

function mix_models(σs, ϕ, ::single_phase)
    @assert length(σs) == 1
    return first(σs)
end
