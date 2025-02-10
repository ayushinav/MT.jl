# We can have a few options on bounds transformation

"""
`sigmoid(m)`:
move to model domain from optimization domain using a sigmoid transformation
"""

function sigmoid(m, bounds) #m, bounds) #m::T1, bounds::Vector{T})::T where {T <: Union{Float32, Float64}, T1}
    σ(x) = 1 / (1 + exp(-x / (bounds[2] - bounds[1])))
    # return 10^(σ(m)*log10(bounds[2]/bounds[1])+ log10(bounds[1]))
    return σ(m) * (bounds[2] - bounds[1]) + bounds[1]
end

function pow_sigmoid(m, bounds) #m::T1, bounds::Vector{T})::T where {T <: Union{Float32, Float64}, T1}
    return 10^ sigmoid(m, bounds)
end

"""
`d_sigmoid(m)`:
gradient for the transformation from optimization domain to model domain. Used for estimating jacobians, but is also useful in analysing sensitivities.
"""
function d_sigmoid(m, bounds) #m::T1, bounds::Vector{T})::T where {T <: Union{Float32, Float64}, T1}
    d_σ(x) = inv((1 + exp(x)) * (1 + exp(-x)))
    # return σ(log10(m))*(log10(bounds[2])- log10(bounds[1]))+ log10(bounds[1])
    return d_σ(m / (bounds[2] - bounds[1]))#*(bounds[2]- bounds[1])
end

function d_pow_sigmoid(m, bounds) #m::T1, bounds::Vector{T})::T where {T <: Union{Float32, Float64}, T1}
    return 10^ sigmoid(m, bounds) * d_sigmoid(m, bounds)
end

"""
`inverse_sigmoid()`: get back to the optimization domain from model domain
"""
function inverse_sigmoid(
        x::T1, bounds::Vector{T})::T where {T <: Union{Float32, Float64}, T1}
    # if x ≈ bounds[2] return 100 end
    # if x ≈ bounds[1] return -100 end

    return (bounds[2] - bounds[1]) * (log(abs(x - bounds[1])) - log(abs(bounds[2] - x)))
end

function inverse_pow_sigmoid(m, bounds) #m::T1, bounds::Vector{T})::T where {T <: Union{Float32, Float64}, T1}
    return inverse_sigmoid(log10(m), bounds)
end

mutable struct transform_utils{T}
    p::Vector{T} # parameters of the transformation function
    tf::Function
    itf::Function
    dtf::Function
end
"""
`transform_utils`: Contains the parameters and functions for transformation from optimization to model domains.
"""
function transform_utils(p::Vector{T}, tf::Function, itf::Function,
        dtf::Function) where {T <: Union{Float32, Float64}}
    return transform_utils{eltype(p)}(
        p, (x) -> tf(x, p), (x) -> itf(x, p), (x) -> dtf(x, p))
end

# should generally be good for most MT inversions
sigmoid_tf = transform_utils([-3.0, 6.0], sigmoid, inverse_sigmoid, d_sigmoid);
pow_tf = transform_utils([], (x) -> 10^x, log10, (x) -> (10^x * log(10)));
log_tf = transform_utils([], log10, (x) -> 10^x, (x) -> inv(x * log(10)));
pow_sigmoid_tf = transform_utils([-3.0, 6.0], pow_sigmoid, inverse_pow_sigmoid, d_pow_sigmoid);
lin_tf = transform_utils([], (x) -> x, (x) -> x, (x) -> 1.0);
