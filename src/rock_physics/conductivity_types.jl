const global boltz_k = 8.617f-5
const global charge_e = 1.602f-19
const global gas_R = 0.008314f0

abstract type AbstractMineralModel end
abstract type AbstractMeltModel end

mutable struct RockphyCond{T} <: AbstractRockphyResponse
    σ::T
end

"""
    SEO3(T)

Laboratory derived relation from naturally occuring Dunite (Olivine)
aggregate, with a small amount of pyroxene pressent to control silica activity.
The SEO3 anhydrous method has a dependency on oxygen fugacity and charge carrier mobility
(as a function of temperature) to calculate the conductivity of Olivine.
The temperature range for the experiment cover 1000 to 1200 degrees C.

## Arguments

  - `T` : Temperature of olivine (in K)

## References

"A new model of olivine electrical conductivity", Geophysical Journal International,
Volume 166, Issue 1, July 2006, Pages 435–437 (https://doi.org/10.1111/j.1365-246X.2006.03041).
"""
mutable struct SEO3{F} <: AbstractMineralModel
    T::F
end

"""
    UHO2014(T, Ch2o_ol)

## Arguments

  - `T` : Temperature of olivine (in K)
  - `Ch2o_ol` : water concentration in olivine (in ppm)

## References

"Toward a unified hydrous olivine electrical conductivity law", [GeoChemistry, Geophysics, Geosystems],
Volume 15, Issue 12 (https://doi.org/10.1002/2014GC005496)
"""
mutable struct UHO2014{F1, F2} <: AbstractMineralModel
    T::F1
    Ch2o_ol::F2
end

# melts
"""
    Ni201(T, Ch2o_m)

Laboratory derived relation from electrical conductivity studies of hydrous and anhydrous basaltic melts.

## Arguments

  - `T` : Temperature of melt (should be greater than 1146.8 K)
  - `Ch2o_m` : water concentration in melt (in ppm)

## References

"Electrical conductivity of hydrous basaltic melts: implications for partial melting in the upper mantle",
Contrib Mineral Petrol 162, 637–650 (2011) (https://doi.org/10.1007/s00410-011-0617-4).
"""
mutable struct Ni2011{F1, F2} <: AbstractMeltModel
    T::F1
    Ch2o_m::F2
    function Ni2011(T, Ch2o_m)
        if T < params_Ni2011.T_corr
            @warn "T (= $T K) should be greater than $(params_Ni2011.T_corr) K for `Ni2011` otherwise erroneous values are obtained"
        end
        return new{typeof(T), typeof(Ch2o_m)}(T, Ch2o_m)
    end
end

# ====== mixing laws =====
struct HS1962_plus end
struct HS1962_minus end
struct HSn_plus end
struct HSn_minus end
struct single_phase end

# == response
