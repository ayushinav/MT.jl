
# """
# Takes in 
# T : temperature (K)
# g : grain size
# Ch2o_ol : water content in olivine (or solid phase (?) ), in ppm
# """
# mutable struct vbr_model{F1, F2, F3} <: AbstractMineralModel
#     T::F1
#     g::F2
#     Ch2o_m::F3
# end

# mutable struct RockphyCond{T} <: AbstractRockphyResponse
#     σ::T
# end
