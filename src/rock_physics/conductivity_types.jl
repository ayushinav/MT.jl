mutable struct RockphyCond{T} <: AbstractRockphyResponse
    σ::T
end

# TODO : use of @concrete types

mutable struct SEO3{F} <: AbstractRockphyModel
    T::F
end

# melts

mutable struct Ni2011{F1, F2}<: AbstractRockphyModel
    T::F1
    Ch2o_m::F2
end

# ====== mixing laws =====
struct HS_plus end
struct HS_minus end
struct single_phase end

# == response




