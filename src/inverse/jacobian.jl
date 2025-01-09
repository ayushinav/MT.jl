"""
`mt_jacobian_cache`: cache `struct` to store the values during the estimation of jacobian.
"""
mutable struct mt_jacobian_cache{T}
    r1::T
    r2::T
end

"""
`mt_jacobian_cache`: cache `struct` to store the values during the estimation of jacobian.
"""
function mt_jacobian_cache(vars) #can change this with zero response
    r1 = forward(MTModel([1.0, 1.0], [1.0]), vars)
    r2 = forward(MTModel([1.0, 1.0], [1.0]), vars)
    return mt_jacobian_cache{typeof(r1)}(r1, r2)
end
"""
`jacobian_mt`: a cache `struct` to store the jacobians, of the same type as the `response` or `response fields` passed.
"""
mutable struct jacobian_mt{T}
    ρₐ::AbstractMatrix{T}
    ϕ::AbstractMatrix{T}
end

"""
`jacobian_mt`: a cache `struct` to store the jacobians, of the same type as the `response` or `response fields` passed.
This is a more general function which will allocate the matrices for only the required fields in both  `model` and `response`.
Useful while testing for different fields of `model` and `response`.
"""
function jacobian_mt(resp_fields, var_eltype, mod, mod_fields)
    empty_mats = [[] for k in resp_fields]
    j = jacobian_mt{var_eltype}(empty_mats...)

    n = sum([length(getfield(mod, k)) for k in mod_fields])

    for k in resp_fields
        setfield!(j, k, zeros(length(ω), n))
    end
    return j
end

"""
`jacobian_mt`: a cache `struct` to store the jacobians, of the same type as the `response` or `response fields` passed.
This is a slightly specific function which will allocate the matrices for only the required fields in `response`.
"""
function jacobian_mt(resp_fields, var_eltype)
    eltypes = [var_eltype for k in resp_fields]
    empty_mats = [[] for k in resp_fields]
    return jacobian_mt{var_eltype}(empty_mats...)
end

# This does not pass a forward function because the cache is mt specific
"""
    jacobian!(J::jacobian_mt,
        m::model, 
        vars::Vector{T},
        mtjc::mt_jacobian_cache;
        model_fields::Vector{Symbol}= [k for k ∈ fieldnames(typeof(m))], 
        response_fields::Vector{Symbol}= [k for k ∈ fieldnames(typeof(J))], 
        ) where T <: Union{Float64, Float32}

overwrites a `jacobian_mt` cache to calculate the jacobian of a `model`. Need to pass `mt_jacobian_cache` for performance.
If no `model_fields` or `response_fields` are passed, all the fields of `model` and the `response` (defined in `jacobian`) will be used.
"""
function jacobian!(J::jacobian_mt,
        m::model,
        vars::Vector{T},
        mtjc::mt_jacobian_cache;
        model_fields::Vector{Symbol}=[k for k in fieldnames(typeof(m))],
        response_fields::Vector{Symbol}=[k for k in fieldnames(typeof(J))]) where {
        T <: Union{Float64, Float32}, model <: AbstractModel}
    nl = 0
    @inbounds for k in model_fields
        ϵ = sqrt(eps(eltype(getfield(m, k))))

        @inbounds for i in eachindex(getfield(m, k))
            getfield(m, k)[i] = getfield(m, k)[i] + ϵ
            MT.forward!(mtjc.r1, m, vars)

            getfield(m, k)[i] = getfield(m, k)[i] - 2ϵ
            forward!(mtjc.r2, m, vars)

            @inbounds for l in response_fields
                getfield(J, l)[:, i] .= getfield(mtjc.r1, l) .- getfield(mtjc.r2, l)
                rmul!(view(getfield(J, l), :, i), inv(2ϵ))
            end
            getfield(m, k)[i] = getfield(m, k)[i] + ϵ
        end
        nl += length(getfield(m, k))
    end
end

# jacobian works fine for 1D MT model, we would have to check for 2D/3D models.
