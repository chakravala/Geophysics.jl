
#   This file is part of Geophysics.jl. It is licensed under the AGPL license
#   Geophysics Copyright (C) 2020 Michael Reed

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

@pure viscosity(G::MoleGas,U::UnitSystem) = viscosity(viscosity(G),U,Metric)
@pure sutherlandviscosity(G::MoleGas,U::UnitSystem) = temperature(sutherlandviscosity(G),U,Metric)

"""
    thermalconductivity(T::Real,G::AbstractMole)

Laminar thermal conductivity `k` of temperature variation (W⋅m⁻¹⋅K⁻¹ or lb⋅s⁻¹⋅°R⁻¹).
""" # 0.5778*778/3600 (English, BTU⋅h?)
@pure thermalconductivity(::MoleGas{M,μ,Tμ,k}) where {M,μ,Tμ,k} = k
@pure sutherlandconductivity(::MoleGas{M,μ,Tμ,k,Tk}) where {M,μ,Tμ,k,Tk} = Tk

@pure thermalconductivity(G::MoleGas,U::UnitSystem) = thermalconductivity(thermalconductivity(G),U,Metric)
@pure sutherlandconductivity(G::MoleGas,U::UnitSystem) = temperature(sutherlandconductivity(G),U,Metric)

gastext(::MoleGas) = "MoleGas{$(molarmass(G)),"

show(io::IO,G::MoleGas) = print(io,gastext(G),"μ=$(viscosity(G)),Tμ=$(sutherlandviscosity(G)),k=$(thermalconductivity(G)),Tk=$(sutherlandconductivity(G))}")

for op ∈ (:viscosity,:thermalconductivity)
    OP = Symbol(:sutherland,op ≠ :viscosity ? :conductivity : op)
    @eval @pure function $op(T::Real,G::MoleGas,U::UnitSystem=Metric)
        Tk = $OP(G,U)
        ((2*$op(G,U))/sqrt(Tk))*(sqrt(T)/(1+Tk/T))
    end
end

@pure function viscond(μ0,Tμ,k0,Tk,T0=288.16)
    T = 2(T0^1.5)
    μ = μ0*sqrt(Tμ)*(T0+Tμ)/T
    k = k0*sqrt(Tk)*(T0+Tk)/T
    return μ,k
end
#=@pure function viscond(G,T=518.69)
    μ = viscosity(T/kelvin,G)*0.020886
    Tμ = sutherlandviscosity(G)/rankine
    k = thermalconductivity(T/kelvin,G)*0.5778*778/3600
    Tk = sutherlandconductivity(G)/rankine
    return μ,Tμ,k,Tk
end=#

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
@pure specificenergy(T::Real,G::AbstractMole,U::UnitSystem=Metric) = heatvolume(T,G,U)*T

"""
    enthalpy(T::Real,G::AbstractMole) = heatpressure(T,G)*T

Specific enthalpy `h` of ideal gas `specificenergy(T,G)+gasconstant(G)*T` (J⋅kg⁻¹ or ft⋅lb⋅slug⁻¹).
"""
@pure specificenthalpy(T::Real,G::AbstractMole,U::UnitSystem=Metric) = heatpressure(T,G,U)*T

"""
    freedom(T::Real,G::AbstractMole) = 2heatvolume(T,G)/gasconstant(G)

Degrees of freedom `f` with storage energy per mole `gasconstant(G)*T/2` (dimensionless).
"""
@pure freedom(T::Real,G::AbstractMole,U::UnitSystem=Metric) = heatvolume(T,G,U)*(2/gasconstant(G,U))

"""
    prandtl(T::Real,G::AbstractMole) = viscosity(T,G)*heatpressure(T,G)/thermalconductivity(T,G)

Prandtl number is the ratio of momentum diffusivity to heat diffusivity (dimensionless).
"""
@pure prandtl(T::Real,G::AbstractMole,U::UnitSystem=Metric) = viscosity(T,G,U)*heatpressure(T,G,U)/thermalconductivity(T,G,U)

"""
    sonicspeed(T::Real,G::AbstractMole) = sqrt(gasconstant(G)*heatratio(T,G)*T)

Speed of sound wave disturbance at isentropic conditions in ideal gas (m⋅s⁻¹ or ft⋅s⁻¹).
"""
@pure sonicspeed(T::Real,G::AbstractMole,U::UnitSystem=Metric) = sqrt((gasconstant(G,U)*heatratio(T,G,U))*T)

export AtomicGas, DiatomicGas, TriatomicGas, SutherlandGas
export wavenumber, frequency, vibration

struct AtomicGas{M,μ,Tμ,k,Tk} <: MoleGas{M,μ,Tμ,k,Tk} end
struct DiatomicGas{M,ν,μ,Tμ,k,Tk} <: MoleGas{M,μ,Tμ,k,Tk} end
struct TriatomicGas{M,ν1,ν2,μ,Tμ,k,Tk} <: MoleGas{M,μ,Tμ,k,Tk} end
struct PentatomicGas{M,μ,Tμ,k,Tk} <: MoleGas{M,μ,Tμ,k,Tk} end
struct SutherlandGas{M,cᵥ,μ,Tμ,k,Tk} <: MoleGas{M,μ,Tμ,k,Tk} end

function AtomicGas(M,μ0,Tμ,k0,Tk,T0=288.16)
    μ,k = viscond(μ0,Tμ,k0,Tk,T0)
    AtomicGas{M,μ,Tμ,k,Tk}()
end
function DiatomicGas(M,ν,μ0,Tμ,k0,Tk,T0=288.16)
    μ,k = viscond(μ0,Tμ,k0,Tk,T0)
    DiatomicGas{M,ν,μ,Tμ,k,Tk}()
end
function TriatomicGas(M,ν1,ν2,μ0,Tμ,k0,Tk,T0=288.16)
    μ,k = viscond(μ0,Tμ,k0,Tk,T0)
    TriatomicGas{M,ν1,ν2,μ,Tμ,k,Tk}()
end
function PentatomicGas(M,μ0,Tμ,k0,Tk,T0=288.16)
    μ,k = viscond(μ0,Tμ,k0,Tk,T0)
    PentatomicGas{M,μ,Tμ,k,Tk}()
end
function SutherlandGas(M,cᵥ,μ0,Tμ,k0,Tk,T0)
    μ,k = viscond(μ0,Tμ,k0,Tk,T0)
    SutherlandGas{M,cᵥ,μ,Tμ,k,Tk}()
end

"""
    wavenumber(G::MoleGas) = frequency(G)/lightspeed(G)

Spectroscopic vibrational wavenumbers of polyatomic molecules if applicable (m⁻¹ or ft⁻¹).
"""
@pure wavenumber(::DiatomicGas{M,ν}) where {M,ν} = ν
@pure wavenumber(::TriatomicGas{M,ν1,ν2}) where {M,ν1,ν2} = ν1,ν2

@pure wavenumber(G::MoleGas,U::UnitSystem) = wavenumber.(wavenumber(G),Ref(U),Ref(Metric))

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
@pure vibration(θT::Float64) = (eθT = exp(θT); θT^2*eθT/(eθT-1)^2)

"""
    heatvolume(T::Real,G::AbstractMole) = translation + rotation + vibration

Specific heat `cᵥ` of ideal gas at constant volume (J⋅kg⁻¹⋅K⁻¹ or ft⋅lb⋅slug⁻¹⋅°R⁻¹).
"""
@pure heatvolume(::Real,G::AtomicGas,U::UnitSystem=Metric) = (3/2)*gasconstant(G,U)
@pure heatvolume(::Real,::SutherlandGas{M,cᵥ}) where {M,cᵥ} = cᵥ
@pure heatvolume(::Real,G::SutherlandGas,U::US) = specificentropy(heatvolume(G),U,Metric)
@pure function heatvolume(T::Real,G::DiatomicGas,U::UnitSystem=Metric)
    R,θ = gasconstant(G,U),vibration(G,U)
    R*((5/2)+vibration(θ/T))
end
@pure function heatvolume(T::Real,G::TriatomicGas,U::UnitSystem=Metric)
    R = gasconstant(G,U)
    θ1,θ2 = vibration(G,U)
    R*((5/2)+vibration(θ1/T)+vibration(θ2/T))
end

gastext(G::AtomicGas) = "AtomicGas{M=$(molarmass(G)),"
gastext(G::DiatomicGas) = "DiatomicGas{M=$(molarmass(G)),ν=$(vibration(G)),"
gastext(G::TriatomicGas) = "TriatomicGas{M=$(molarmass(G)),ν₁=$(vibration(G)[1]),ν₂=$(vibration(G)[2]),"
gastext(G::SutherlandGas) = "Gas{M=$(molarmass(G)),cᵥ=$(heatvolume(G)),cₚ=$(heatpressure(G)),"

struct Mixture{M,N,C} <: AbstractMole{M}
    f::Values{N,Float64}
end

Mixture{N,C}(f) where {N,C} = Mixture{f⋅relativemass.(C),N,C}(f)
Mixture{C}(f) where {C,U} = Mixture{length(C),C}(f)

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

