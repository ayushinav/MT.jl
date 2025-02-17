"""
`jacobian_cache`: cache `struct` to store the values during the estimation of jacobian.
"""
mutable struct jacobian_cache{T1, T2, T3}
    r1::T1
    r2::T2
    j::T3
end

function jacobian_cache(resp_fields, resp, mod, mod_fields)
    arr = []
    for k in resp_fields
        for j in mod_fields
            push!(arr, zeros(length(getfield(resp, k)), length(getfield(mod, j))))
        end
    end

    resp_type = (typeof(resp)).name.wrapper{[AbstractMatrix for k in resp_fields]...}

    return jacobian_cache(zero(resp), zero(resp), resp_type(arr...))
end

"""
    jacobian!(m::model, 
        vars::Vector{T},
        jc::jacobian_cache;
        model_fields::Vector{Symbol}= [k for k ∈ fieldnames(typeof(m))], 
        response_fields::Vector{Symbol}= [k for k ∈ fieldnames(typeof(J))], 
        response_trans_utils::NamedTuple=(; ρₐ=lin_tf, ϕ=lin_tf)
        ) where T <: Union{Float64, Float32}

overwrites a `jacobian_cache` cache to calculate the jacobian of a `model`.
If no `model_fields` or `response_fields` are passed, all the fields of `model` and the `response` (defined in `jacobian`) will be used.
"""
function jacobian!(m::model,
        vars::Vector{T},
        jc::jacobian_cache;
        model_fields::Vector{Symbol}=[k for k in fieldnames(typeof(m))],
        response_fields::Vector{Symbol}=[k for k in fieldnames(typeof(J))],
        response_trans_utils::NamedTuple=(; ρₐ=lin_tf, ϕ=lin_tf)) where {
        T <: Union{Float64, Float32}, model <: AbstractModel}
    nl = 0
    @inbounds for k in model_fields
        ϵ = sqrt(eps(eltype(getfield(m, k))))

        @inbounds for i in eachindex(getfield(m, k))
            getfield(m, k)[i] = getfield(m, k)[i] + ϵ
            forward!(jc.r1, m, vars; response_trans_utils=response_trans_utils)

            getfield(m, k)[i] = getfield(m, k)[i] - 2ϵ
            forward!(jc.r2, m, vars; response_trans_utils=response_trans_utils)

            @inbounds for l in response_fields
                view(getfield(jc.j, l), :, i) .= getfield(jc.r1, l) .- getfield(jc.r2, l)
                rmul!(view(getfield(jc.j, l), :, i), inv(2ϵ))
            end
            getfield(m, k)[i] = getfield(m, k)[i] + ϵ
        end
        nl += length(getfield(m, k))
    end
end

# jacobian works fine for 1D MT model, we would have to check for 2D/3D models.
