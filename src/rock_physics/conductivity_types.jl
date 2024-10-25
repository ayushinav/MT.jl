mutable struct RockphyCond{T} <: AbstractRockphyResponse
    m::T
end

# TODO : use of @concrete types

mutable struct SEO3{F} <: AbstractRockphyModel
    T::F
end


mutable struct model2{F1, F2} <: AbstractRockphyModel
    T::F
    ϕ::F2
end

mutable struct model2{F1, F2} <: AbstractRockphyModel
    T::F
    ϕ::F2
end


# melts

mutable struct Ni2011{F1, F2}<: AbstractRockphyModel
    T::F1
    Ch2o_m::F2
end


# ====== mixing laws =====

mutable struct HS_plus{T1 <: AbstractRockphyModel, T2} <: AbstractRockphyModel
    models::Vector{T1}
    ϕ::Vector{T2}

    # test that ∑ϕ = 1
end


# == response




