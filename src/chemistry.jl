
#   This file is part of Geophysics.jl
#   It is licensed under the AGPL license
#   Geophysics Copyright (C) 2020 Michael Reed
#       _           _                         _
#      | |         | |                       | |
#   ___| |__   __ _| | ___ __ __ ___   ____ _| | __ _
#  / __| '_ \ / _` | |/ / '__/ _` \ \ / / _` | |/ _` |
# | (__| | | | (_| |   <| | | (_| |\ V / (_| | | (_| |
#  \___|_| |_|\__,_|_|\_\_|  \__,_| \_/ \__,_|_|\__,_|
#
#   https://github.com/chakravala
#   https://crucialflow.com

export molarmass, units, molecularmass, gasconstant, AbstractMole, MoleGas
export sutherlandviscosity, sutherlandconductivity

# fluid models

"""
    abstract type AbstractMole{M} end

Chemical susbtance specified by molarmass `M`.
Induces dervived values of `molarmass`, `units`, `molecularmass`, `gasconstant`,`viscosity`, `thermalconductivity`, `heatvolume`, `heatpressure`, `heatratio`, `specificenergy`, `specificenthalpy`, `freedom`, `prandtl`, `sonicspeed`, and constants associated with the `UnitSystem` values.
"""
abstract type AbstractMole{M} end

"""
    relativemass(x) = molecularmass(x)/atomicmass(units(x))

Relative atomic mass is ratio of average mass per (atom) to `atomicmass` (dimensionless).
"""
@pure relativemass(::AbstractMole{M}) where M = M

"""
    molarmass(x) = molarmass(units(x))*relativemass(x)

Mass of 1 mole of a chemical `molecularmass(x)*avogadro(x)` (kg⋅mol⁻¹ or slug⋅slug-mol⁻¹).
"""
@pure molarmass(M::AbstractMole,U::UnitSystem=Metric) = molarmass(U)*relativemass(M)

"""
    molecularmass(x) = molarmass(x)/avogadro(x)

Mass of 1 molecule of a chemical (kg or slugs).
"""
@pure molecularmass(M::AbstractMole,U::UnitSystem=Metric) = molarmass(M,U)/avogadro(U)

"""
    gasconstant(x) = universal(x)/molarmass(x)

Specific gas constant of chemical (J⋅kg⁻¹⋅K⁻¹ or ft⋅lb⋅slug⁻¹⋅°R⁻¹).
"""
@pure gasconstant(M::AbstractMole,U::UnitSystem=Metric) = universal(U)/molarmass(M,U)

for op ∈ Constants
    op ≠ :molarmass && @eval begin; export $op
        @pure $op(M::AbstractMole,U::UnitSystem=Metric) = $op(U)
    end
end

"""
    abstract type MoleGas{M,μ,Tμ,k,Tk} <: AbstractMole{M} end

Molecular gas susbtance specified by molarmass `M` and Sutherland constants `μ,Tμ,k,Tk`.
Induces additional dervived values if applicable `wavenumber`, `wavelength`, `frequency`, `vibration`, and inherited constants associated with the `AbstractMole` values.
"""
abstract type MoleGas{M,μ,Tμ,k,Tk} <: AbstractMole{M} end

"""
    viscosity(T::Real,G::AbstractMole)

Laminar dynamic vicsosity `μ` is stress to normal acceleration ratio (Pa⋅s or slug⋅ft⁻¹⋅s⁻¹).
"""
@pure viscosity(::MoleGas{M,μ}) where {M,μ} = μ
@pure sutherlandviscosity(::MoleGas{M,μ,Tμ}) where {M,μ,Tμ} = Tμ

@pure viscosity(G::MoleGas,U::UnitSystem) = viscosity(G)*viscosity(Metric,U)
@pure sutherlandviscosity(G::MoleGas,U::UnitSystem) = sutherlandviscosity(G)*temperature(Metric,U)

"""
    thermalconductivity(T::Real,G::AbstractMole)

Laminar thermal conductivity `k` of temperature variation (W⋅m⁻¹⋅K⁻¹ or lb⋅s⁻¹⋅°R⁻¹).
""" # 0.5778*778/3600 (English, BTU⋅h?)
@pure thermalconductivity(::MoleGas{M,μ,Tμ,k}) where {M,μ,Tμ,k} = k
@pure sutherlandconductivity(::MoleGas{M,μ,Tμ,k,Tk}) where {M,μ,Tμ,k,Tk} = Tk

@pure thermalconductivity(G::MoleGas,U::UnitSystem) = thermalconductivity(G)*thermalconductivity(Metric,U)
@pure sutherlandconductivity(G::MoleGas,U::UnitSystem) = sutherlandconductivity(G)*temperature(Metric,U)

gastext(::MoleGas) = "MoleGas{$(molarmass(G)),"

show(io::IO,G::MoleGas) = print(io,gastext(G),"μ=$(viscosity(G)),Tμ=$(sutherlandviscosity(G)),k=$(thermalconductivity(G)),Tk=$(sutherlandconductivity(G))}")

for op ∈ (:viscosity,:thermalconductivity)
    OP = Symbol(:sutherland,op ≠ :viscosity ? :conductivity : op)
    @eval @pure function $op(T::Real,G::MoleGas,U::UnitSystem=Metric)
        Tk = $OP(G,U)
        ((2*$op(G,U))/sqrt(normal(Tk)))*(sqrt(T)/(1+normal(Tk)/T))
    end
end

@pure function viscond(μ0,Tμ,k0,Tk,T0=288.16,U=Metric)
    T1 = 2(T0^1.5)
    μ = μ0*sqrt(Tμ)*(T0+Tμ)/T1
    k = k0*sqrt(Tk)*(T0+Tk)/T1
    return Quantity(F*T/L/L,U,μ),Quantity(F/T/Θ,U,k)
end

for heat ∈ (:heatratio,:heatvolume,:heatpressure)
    @eval @pure $heat(M::AbstractMole,U::US=Metric) = $heat(temperature(288.16,U,Metric),M,U)
end

"""
    heatpressure(T::Real,G::AbstractMole) = heatvolume(T,G)+gasconstant(G)

Specific heat `cₚ` of ideal gas at constant pressure (J⋅kg⁻¹⋅K⁻¹ or ft⋅lb⋅slug⁻¹⋅°R⁻¹).
"""
@pure heatpressure(T::Real,G::MoleGas,U::UnitSystem=Metric) = heatvolume(T,G,U)+gasconstant(G,U)

"""
    heatratio(T::Real,G::AbstractMole) = heatpressure(T,G)/heatvolume(T,G)

Specific heat ratio `γ` at constant pressure to constant volume of ideal gas (dimensionless).
"""
@pure heatratio(T::Real,G::MoleGas,U::UnitSystem=Metric) = (gasconstant(G,U)/heatvolume(T,G,U))+1

"""
    specificenergy(T::Real,G::AbstractMole) = heatvolume(T,G)*T

Specific energy `e` of ideal gas `specificenthalpy(T,G)-gasconstant(G)*T` (J⋅kg⁻¹ or ft⋅lb⋅slug⁻¹).
"""
@pure specificenergy(T::Real,G::AbstractMole,U::UnitSystem=Metric) = heatvolume(T,G,U)*Quantity(Θ,U,T)

"""
    enthalpy(T::Real,G::AbstractMole) = heatpressure(T,G)*T

Specific enthalpy `h` of ideal gas `specificenergy(T,G)+gasconstant(G)*T` (J⋅kg⁻¹ or ft⋅lb⋅slug⁻¹).
"""
@pure specificenthalpy(T::Real,G::AbstractMole,U::UnitSystem=Metric) = heatpressure(T,G,U)*Quantity(Θ,U,T)

"""
    freedom(T::Real,G::AbstractMole) = 2heatvolume(T,G)/gasconstant(G)

Degrees of freedom `f` with storage energy per mole `gasconstant(G)*T/2` (dimensionless).
"""
@pure freedom(T::Real,G::AbstractMole,U::UnitSystem=Metric) = heatvolume(T,G,U)*(2/gasconstant(G,U))

"""
    prandtl(T::Real,G::AbstractMole) = viscosity(T,G)*heatpressure(T,G)/thermalconductivity(T,G)

Prandtl number is the ratio of momentum diffusivity to heat diffusivity (dimensionless).
"""
@pure prandtl(T::Real,G::AbstractMole,U::UnitSystem=Metric) = gravity(U)*viscosity(T,G,U)*heatpressure(T,G,U)/thermalconductivity(T,G,U)

"""
    sonicspeed(T::Real,G::AbstractMole) = sqrt(gasconstant(G)*heatratio(T,G)*T)

Speed of sound wave disturbance at isentropic conditions in ideal gas (m⋅s⁻¹ or ft⋅s⁻¹).
"""
@pure sonicspeed(T::Real,G::AbstractMole,U::UnitSystem=Metric) = sqrt((gasconstant(G,U)*gravity(U)*heatratio(T,G,U))*Quantity(Θ,U,T))

export AtomicGas, DiatomicGas, TriatomicGas, SutherlandGas
export wavenumber, frequency, vibration

struct AtomicGas{M,μ,Tμ,k,Tk} <: MoleGas{M,μ,Tμ,k,Tk} end
struct DiatomicGas{M,ν,μ,Tμ,k,Tk} <: MoleGas{M,μ,Tμ,k,Tk} end
struct TriatomicGas{M,ν1,ν2,μ,Tμ,k,Tk} <: MoleGas{M,μ,Tμ,k,Tk} end
struct PentatomicGas{M,μ,Tμ,k,Tk} <: MoleGas{M,μ,Tμ,k,Tk} end
struct SutherlandGas{M,cᵥ,μ,Tμ,k,Tk} <: MoleGas{M,μ,Tμ,k,Tk} end

function AtomicGas(M,μ0,Tμ,k0,Tk,T0=288.16,U=Metric)
    μ,k = viscond(μ0,Tμ,k0,Tk,T0,U)
    AtomicGas{M,μ,Quantity(Θ,U,Tμ),k,Quantity(Θ,U,Tk)}()
end
function DiatomicGas(M,ν,μ0,Tμ,k0,Tk,T0=288.16,U=Metric)
    μ,k = viscond(μ0,Tμ,k0,Tk,T0,U)
    DiatomicGas{M,Quantity(inv(L),U,ν),μ,Quantity(Θ,U,Tμ),k,Quantity(Θ,U,Tk)}()
end
function TriatomicGas(M,ν1,ν2,μ0,Tμ,k0,Tk,T0=288.16,U=Metric)
    μ,k = viscond(μ0,Tμ,k0,Tk,T0,U)
    TriatomicGas{M,Quantity(inv(L),U,ν1),Quantity(inv(L),U,ν2),μ,Quantity(Θ,U,Tμ),k,Quantity(Θ,U,Tk)}()
end
function PentatomicGas(M,μ0,Tμ,k0,Tk,T0=288.16,U=Metric)
    μ,k = viscond(μ0,Tμ,k0,Tk,T0,U)
    PentatomicGas{M,μ,Quantity(Θ,U,Tμ),k,Quantity(Θ,U,Tk)}()
end
function SutherlandGas(M,cᵥ,μ0,Tμ,k0,Tk,T0=288.16,U=Metric)
    μ,k = viscond(μ0,Tμ,k0,Tk,T0,U)
    SutherlandGas{M,Quantity(F*L/M/Θ,U,cᵥ),μ,Quantity(Θ,U,Tμ),k,Quantity(Θ,U,Tk)}()
end

"""
    wavenumber(G::MoleGas) = frequency(G)/lightspeed(G)

Spectroscopic vibrational wavenumbers of polyatomic molecules if applicable (m⁻¹ or ft⁻¹).
"""
@pure wavenumber(::DiatomicGas{M,ν}) where {M,ν} = ν
@pure wavenumber(::TriatomicGas{M,ν1,ν2}) where {M,ν1,ν2} = ν1,ν2

@pure wavenumber(G::MoleGas,U::UnitSystem) = wavenumber(G).*Ref(wavenumber(Metric,U))

"""
    wavelength(G::MoleGas) = 1/wavenumber(G)

Spectroscopic vibrational wavelength `λ` of polyatomic molecules if applicable (m or ft).
"""
@pure wavelength(G::MoleGas,U::UnitSystem=Metric) = inv.(wavenumber(G,U))

"""
    frequency(G::MoleGas) = wavenumber(G)*lightspeed(G)

Spectroscopic vibrational frequencies `ν` of polyatomic molecules if applicable (Hz or s⁻¹).
"""
@pure frequency(G::MoleGas,U::UnitSystem=Metric) = wavenumber(G,U).*lightspeed(U)

"""
    vibration(G::MoleGas) = frequency(G)*planck(G)/boltzmann(G)

Vibrational temperature `θ` of polyatomic molecules if applicable (K or °R).
"""
@pure vibration(G::MoleGas,U::US=Metric) = frequency(G,U).*(planck(U)/boltzmann(U)/1.2)
@pure vibration(θT::T) where T<:AbstractFloat = (eθT = exp(θT); θT^2*eθT/(eθT-1)^2)
if usingSimilitude
@pure vibration(θT::Quantity{D,U}) where {D,U} = Quantity{D,U}(vibration(θT.v))
end

"""
    heatvolume(T::Real,G::AbstractMole) = translation + rotation + vibration

Specific heat `cᵥ` of ideal gas at constant volume (J⋅kg⁻¹⋅K⁻¹ or ft⋅lb⋅slug⁻¹⋅°R⁻¹).
"""
@pure heatvolume(::Real,G::AtomicGas,U::UnitSystem=Metric) = (3/2)*gasconstant(G,U)
@pure heatvolume(::Real,::SutherlandGas{M,cᵥ}) where {M,cᵥ} = cᵥ
@pure heatvolume(::Real,G::SutherlandGas,U::US) = specificentropy(heatvolume(G),U,Metric)
@pure function heatvolume(T::Real,G::DiatomicGas,U::UnitSystem=Metric)
    R,θ = gasconstant(G,U),vibration(G,U)
    R*((5/2)+vibration(normal(θ)/T))
end
@pure function heatvolume(T::Real,G::TriatomicGas,U::UnitSystem=Metric)
    R = gasconstant(G,U)
    θ1,θ2 = vibration(G,U)
    R*((5/2)+vibration(normal(θ1)/T)+vibration(normal(θ2)/T))
end

gastext(G::AtomicGas) = "AtomicGas{M=$(molarmass(G)),"
gastext(G::DiatomicGas) = "DiatomicGas{M=$(molarmass(G)),ν=$(vibration(G)),"
gastext(G::TriatomicGas) = "TriatomicGas{M=$(molarmass(G)),ν₁=$(vibration(G)[1]),ν₂=$(vibration(G)[2]),"
gastext(G::SutherlandGas) = "Gas{M=$(molarmass(G)),cᵥ=$(heatvolume(G)),cₚ=$(heatpressure(G)),"

struct Mixture{M,N,C} <: AbstractMole{M}
    f::Values{N,Float64}
end

Mixture{N,C}(f) where {N,C} = Mixture{f⋅relativemass.(C),N,C}(f)
Mixture{C}(f) where C = Mixture{length(C),C}(f)

@pure chemical(::Mixture{M,N,C}) where {M,N,C} = C
@pure molecules(::Mixture{M,N,C}) where {M,N,C} = Values(C)
@pure fractions(M::Mixture) = M.f
@pure Base.length(::Mixture{M,N}) where {M,N} = N
@pure Base.getindex(M::Mixture,i::Int) = chemical(M)[i]
@pure heatratio(T::Real,M::Mixture,U::UnitSystem=Metric) = heatpressure(T,M,U)/heatvolume(T,M,U)

for op ∈ (:viscosity,:thermalconductivity,:heatvolume,:heatpressure)
    @eval $op(T::Real,M::Mixture,U::UnitSystem=Metric) = fractions(M)⋅Values($op.(T,chemical(M),Ref(U)))
end
for op ∈ (:viscosity,:thermalconductivity,:sutherlandviscosity,:sutherlandconductivity)
    @eval $op(M::Mixture,U::UnitSystem=Metric) = fractions(M)⋅Values($op.(chemical(M),Ref(U)))
end

Base.:*(f::Float64,G::AbstractMole{M}) where M = Mixture{M,1,(G,)}(Values(f))

function Base.:+(m::Mixture...)
    N = sum(length.(m))
    Mixture{N,Tuple(vcat(molecules.(m)...))}(Values{N}(vcat(fractions.(m)...)))
end

# thermodynamic states

"""
    FluidState{f}

Thermodynamic state of fluid `f` at temperature `T` and pressure `P`.
Induces derived values `fluid`, `temperature`, `pressure`, `density`, `specificvolume`, `kinematic`, `heatcapacity`, `thermaldiffusivity`, `elasticity`, `specificimpedance`, `intensity`, and values associated with `f::AbstractMole` derivations.
"""
struct FluidState{f,u}
    T::Float64
    P::Float64
end

for Gas ∈ (:SutherlandGas,:Mixture,:AtomicGas,:DiatomicGas,:TriatomicGas)
    @eval (G::$Gas)(T=288.15,P=atm,U=Metric) = FluidState{G,U}(T,P)
end

@pure units(::FluidState{f,u}) where {f,u} = u
(U::UnitSystem)(F::FluidState{G,S}) where {G,S} = FluidState{G,U}(temperature(F,U),pressure(F,U))

"""
    fluid(x)

Specification of an `AbstractMole` fluid chemical instance.
"""
@pure fluid(::FluidState{f}) where f = f

"""
    temperature(F::FluidState)

Absolute thermodynamic temperature scale `T` at `F` (K or °R).
"""
@pure temperature(F::FluidState) = F.T
@pure temperature(F::FluidState,U::UnitSystem) = temperature(temperature(F),U,units(F))

"""
    pressure(F::FluidState)

Absolute force per unit area `P` exerted by `F` (Pa or lb⋅ft⁻²).
"""
@pure pressure(F::FluidState) = F.P
@pure pressure(F::FluidState,U::UnitSystem) = pressure(pressure(F),U,units(F))

for op ∈ Properties
    @eval @pure $op(F::FluidState,U::US=units(F)) = $op(fluid(F),U)
end
for op ∈ Intrinsic
    @eval @pure $op(F::FluidState,U::US=units(F)) = $op(temperature(F,U),fluid(F),U)
end

@doc """
    viscosity(F::FluidState) = viscosity(temperature(F),fluid(F))

Laminar dynamic vicsosity `μ` at the temperature of `F` (Pa⋅s or lb⋅s¹⋅ft⁻²).
""" viscosity(F::FluidState)

@doc """
    thermalconductivity(F::FluidState) = conductivity(temperature(F),fluid(F))

Laminar thermal conductivity `k` at the temperature `F` (W⋅m⁻¹⋅K⁻¹ or lb⋅s⁻¹⋅°R⁻¹).
""" thermalconductivity(F::FluidState)

@doc """
    heatvolume(F::FluidState) = heatvolume(temperature(F),fluid(F))

Specific heat `cᵥ` at the temperature of `F` (J⋅kg⁻¹⋅K⁻¹ or ft⋅lb⋅slug⁻¹⋅°R⁻¹).
""" heatvolume(F::FluidState)

@doc """
    heatpressure(F::FluidState) = heatpressure(temperature(F),fluid(F))

Specific heat `cₚ` at the temperature of `F` (J⋅kg⁻¹⋅K⁻¹ or ft⋅lb⋅slug⁻¹⋅°R⁻¹).
""" heatpressure(F::FluidState)

@doc """
    heatratio(F::FluidState) = heatratio(temperature(F),fluid(F))

Specific heat ratio `γ` at the temperature of `F` (dimensionless).
""" heatratio(F::FluidState)

@doc """
    specificenergy(F::FluidState) = specificenergy(temperature(F),fluid(F))

Specific energy `e` at the temperature of `F` (J⋅kg⁻¹ or ft⋅lb⋅slug⁻¹).
""" specificenergy(F::FluidState)

@doc """
    specificenthalpy(F::FluidState) = specificenthalpy(temperature(F),fluid(F))

Specific enthalpy `h` at the temperature of `F` (J⋅kg⁻¹ or ft⋅lb⋅slug⁻¹).
""" specificenthalpy(F::FluidState)

@doc """
    freedom(F::FluidState) = freedom(temperature(F),fluid(F))

Degrees of freedom `f` at the temperature of `F` (dimensionless).
""" freedom(F::FluidState)

@doc """
    prandtl(F::FluidState) = prandtl(temperature(F),fluid(F))

Prandtl number ratio at the temperature of `F` (dimensionless).
""" prandtl(F::FluidState)

@doc """
    sonicspeed(F::FluidState) = sonicspeed(temperature(F),fluid(F))

Speed of sound wave disturbance at the temperature of `F` (m⋅s⁻¹ or ft⋅s⁻¹).
""" sonicspeed(F::FluidState)

"""
    density(F::FluidState) = pressure(F)/temperature(F)/gasconstant(F)

Inertial mass per volume `ρ` of at a pressure and temperature (kg⋅m⁻³ or slugs⋅ft⁻³).
"""
@pure density(F::FluidState,U::US=units(F)) = (pressure(F,U)/temperature(F,U))/gasconstant(F,U)

"""
    specificvolume(F::FluidState) = 1/density(F)

Specific volume per mass `v` at a pressure and temperature (m³⋅kg⁻¹ or ft³⋅slug⁻¹).
"""
@pure specificvolume(F::FluidState,U::US=units(F)) = inv(density(F,U))

"""
    kinematic(F::FluidState) = viscosity(F)/density(F)

Kinematic viscosity ratio `ν` at a pressure and temperature (m²⋅s⁻¹ or ft²⋅s⁻¹).
"""
@pure kinematic(F::FluidState,U::US=units(F)) = viscosity(F,U)/density(F,U)

"""
    heatcapacity(F::FluidState) = heatpressure(F)*density(F)

Specific heat per mass at a pressure and temperature (J⋅m⁻³⋅K⁻¹ or lb⋅ft⁻²⋅°R⁻¹).
"""
@pure heatcapacity(F::FluidState,U::US=units(F)) = heatpressure(F,U)*density(F,U)

"""
    thermaldiffusivity(F::FluidState) = thermalconductivity(F)/heatcapacity(F)

Thermal diffusivity `α` at a pressure and temperature (m²⋅s⁻¹ or ft²⋅s⁻¹).
"""
@pure thermaldiffusivity(F::FluidState,U::US=units(F)) = thermalconductivity(F,U)/heatcapacity(F,U)

"""
    elasticity(F::FluidState) = heatratio(F)*pressure(F)

Bulk modulus of elasticity `B` at a pressure and temperature (Pa or slug⋅ft⁻¹⋅s⁻²).
"""
@pure elasticity(F::FluidState,U::US=units(F)) = heatratio(F,U)*pressure(F,U)

"""
    specificimpedance(F::FluidState) = density(F)*sonicspeed(F)

Specific acoustic resistance at a pressure and temperature (kg⋅m⁻³⋅s⁻¹ or slug⋅ft⁻³⋅s⁻¹).
"""
@pure specificimpedance(F::FluidState,U::US=units(F)) = density(F,U)*sonicspeed(F,U)

"""
    intensity(F::FluidState) = pressure(F)^2/impedance(F)

Instantaneous acoustic intensity `I` at a pressure and temperature (W⋅m⁻² or slug⋅s⁻³).
""" # irradiance
@pure intensity(F::FluidState,U::US=units(F)) = pressure(F,U)^2/impedance(F,U)
