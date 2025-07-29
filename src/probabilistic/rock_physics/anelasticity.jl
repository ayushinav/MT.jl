struct RockPhyanelasticDistribution{
    T1 <: Union{Distribution, AbstractArray}, T2 <: Union{Distribution, AbstractArray},
    T3 <: Union{Distribution, AbstractArray}, T4 <: Union{Distribution, AbstractArray},
    T5 <: Union{Distribution, AbstractArray}, T6 <: Union{Distribution, AbstractArray}} <:
       AbstractRockphyResponseDistribution
    J1::T1
    J2::T2
    Qinv::T3
    M::T4
    V::T5
    Vave::T6
end

"""
    andrade_pspDistribution(T, P, dg, σ, ϕ, ρ, f)

Model Distribution for `andrade_psp`[@ref].

## Arguments

    - `T` : Temperature of the rock (K)
    - `P` : Pressure (GPa)
    - `dg`: Grain size (μm)
    - `σ` : Shear stress (GPa)
    - `ϕ` : Porosity
    - `ρ` : Density (kg/m³)
    - `f` : frequency

## Usage

Refer to the documentation for usage examples.
"""
mutable struct andrade_pspDistribution{
    T1 <: Union{Distribution, AbstractArray}, T2 <: Union{Distribution, AbstractArray},
    T3 <: Union{Distribution, AbstractArray}, T4 <: Union{Distribution, AbstractArray},
    T5 <: Union{Distribution, AbstractArray}, T6 <: Union{Distribution, AbstractArray},
    T7 <: Union{Distribution, AbstractArray}} <: AbstractRockphyModelDistribution
    T::T1
    P::T2
    dg::T3
    σ::T4
    ϕ::T5
    ρ::T6
    f::T7
end

"""
    eburgers_pspDistribution(T, P, dg, σ, ϕ, ρ, f)

Model Distribution for `eburgers_psp`[@ref].

## Arguments

    - `T` : Temperature of the rock (K)
    - `P` : Pressure (GPa)
    - `dg`: Grain size (μm)
    - `σ` : Shear stress (GPa)
    - `ϕ` : Porosity
    - `ρ` : Density (kg/m³)
    - `f` : frequency

## Usage

Refer to the documentation for usage examples.
"""
mutable struct eburgers_pspDistribution{
    T1 <: Union{Distribution, AbstractArray}, T2 <: Union{Distribution, AbstractArray},
    T3 <: Union{Distribution, AbstractArray}, T4 <: Union{Distribution, AbstractArray},
    T5 <: Union{Distribution, AbstractArray}, T6 <: Union{Distribution, AbstractArray},
    T7 <: Union{Distribution, AbstractArray}, T8 <: Union{Distribution, AbstractArray},
    T9 <: Union{Distribution, AbstractArray}} <: AbstractRockphyModelDistribution
    T::T1
    P::T2
    dg::T3
    σ::T4
    ϕ::T5
    ρ::T6
    Ch2o_ol::T7
    T_solidus::T8
    f::T9
end

"""
    premelt_anelasticDistribution(T, P, dg, σ, ϕ, ρ, f)

Model Distribution for `premelt_anelastic`[@ref].

## Arguments

    - `T` : Temperature of the rock (K)
    - `P` : Pressure (GPa)
    - `dg`: Grain size (μm)
    - `σ` : Shear stress (GPa)
    - `ϕ` : Porosity
    - `ρ` : Density (kg/m³)
    - `f` : frequency

## Usage

Refer to the documentation for usage examples.
"""
mutable struct premelt_anelasticDistribution{
    T1 <: Union{Distribution, AbstractArray}, T2 <: Union{Distribution, AbstractArray},
    T3 <: Union{Distribution, AbstractArray}, T4 <: Union{Distribution, AbstractArray},
    T5 <: Union{Distribution, AbstractArray}, T6 <: Union{Distribution, AbstractArray},
    T7 <: Union{Distribution, AbstractArray}, T8 <: Union{Distribution, AbstractArray},
    T9 <: Union{Distribution, AbstractArray}} <: AbstractRockphyModelDistribution
    T::T1
    P::T2
    dg::T3
    σ::T4
    ϕ::T5
    ρ::T6
    Ch2o_ol::T7
    T_solidus::T8
    f::T9
end

"""
    xfit_mxwDistribution(T, P, dg, σ, ϕ, ρ, f)

Model Distribution for `xfit_mxw`[@ref].

## Arguments

    - `T` : Temperature of the rock (K)
    - `P` : Pressure (GPa)
    - `dg`: Grain size (μm)
    - `σ` : Shear stress (GPa)
    - `ϕ` : Porosity
    - `ρ` : Density (kg/m³)
    - `f` : frequency

## Usage

Refer to the documentation for usage examples.
"""
mutable struct xfit_mxwDistribution{
    T1 <: Union{Distribution, AbstractArray}, T2 <: Union{Distribution, AbstractArray},
    T3 <: Union{Distribution, AbstractArray}, T4 <: Union{Distribution, AbstractArray},
    T5 <: Union{Distribution, AbstractArray}, T6 <: Union{Distribution, AbstractArray},
    T7 <: Union{Distribution, AbstractArray}, T8 <: Union{Distribution, AbstractArray},
    T9 <: Union{Distribution, AbstractArray}} <: AbstractRockphyModelDistribution
    T::T1
    P::T2
    dg::T3
    σ::T4
    ϕ::T5
    ρ::T6
    Ch2o_ol::T7
    T_solidus::T8
    f::T9
end

"""
    andrade_analyticalDistribution(T, P, dg, σ, ϕ, ρ, f)

Model Distribution for `xfit_mxw`[@ref].

## Arguments

    - `T` : Temperature of the rock (K)
    - `P` : Pressure (GPa)
    - `dg`: Grain size (μm)
    - `σ` : Shear stress (GPa)
    - `ϕ` : Porosity
    - `ρ` : Density (kg/m³)
    - `Ch2o_ol` : water concentration in olivine (in ppm), only used when using `HK2003` for viscosity calculations
    - `T_solidus` : Solidus temperature (K), only used when using `xfit_premelt` for viscosity calculations
    - `f` : frequency

## Usage

Refer to the documentation for usage examples.
"""
mutable struct andrade_analyticalDistribution{
    T1 <: Union{Distribution, AbstractArray}, T2 <: Union{Distribution, AbstractArray},
    T3 <: Union{Distribution, AbstractArray}, T4 <: Union{Distribution, AbstractArray},
    T5 <: Union{Distribution, AbstractArray}, T6 <: Union{Distribution, AbstractArray},
    T7 <: Union{Distribution, AbstractArray}, T8 <: Union{Distribution, AbstractArray},
    T9 <: Union{Distribution, AbstractArray}} <: AbstractRockphyModelDistribution
    T::T1
    P::T2
    dg::T3
    σ::T4
    ϕ::T5
    ρ::T6
    Ch2o_ol::T7
    T_solidus::T8
    f::T9
end
