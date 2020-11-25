
#   This file is part of Geophysics.jl. It is licensed under the AGPL license
#   Geophysics Copyright (C) 2020 Michael Reed

export molarmass, units, molecularmass, gasconstant, AbstractMole, MoleGas
export sutherlandviscosity, sutherlandconductivity

# fluid models

"""
    abstract type AbstractMole{M,U} end

Chemical susbtance specified by molarmass `M` with unit system `U`.
Induces dervived values of `molarmass`, `units`, `molecularmass`, `gasconstant`,`viscosity`, `conductivity`, `heatvolume`, `heatpressure`, `heatratio`, `energy`, `enthalpy`, `freedom`, `prandtl`, `sonicspeed`, and constants associated with the `UnitSystem` values.
"""
abstract type AbstractMole{M,U} end

"""
    molarmass(x) = molecularmass(x)*avogadro(x)

Mass of 1 mole of a chemical (kg⋅mol⁻¹ or slug⋅slug-mol⁻¹).
"""
@pure molarmass(::AbstractMole{M}) where M = M

"""
    units(x)::UnitSystem

Representation of unit system in terms of foundational constants (`Metric` or `English`).
"""
@pure units(::AbstractMole{M,U}) where {M,U} = U

"""
    molecularmass(x) = molarmass(x)/avogadro(x)

Mass of 1 molecule of a chemical (kg or slugs).
"""
@pure molecularmass(F::AbstractMole) = molarmass(F)/avogadro(F)

"""
    gasconstant(x) = universal(x)/molarmass(x)

Specific gas constant of chemical (J⋅kg⁻¹⋅K⁻¹ or ft⋅lb⋅slug⁻¹⋅°R⁻¹).
"""
@pure gasconstant(F::AbstractMole) = universal(F)/molarmass(F)

for op ∈ Constants
    @eval begin; export $op
        @pure $op(F::AbstractMole) = $op(units(F))
    end
end

"""
    abstract type MoleGas{M,μ,Tμ,k,Tk,U} <: AbstractMole{M,U} end

Molecular gas susbtance specified by molarmass `M` and Sutherland constants `μ,Tμ,k,Tk` with unit system `U`.
Induces additional dervived values if applicable `wavenumber`, `wavelength`, `frequency`, `vibration`, and inherited constants associated with the `AbstractMole` values.
"""
abstract type MoleGas{M,μ,Tμ,k,Tk,U} <: AbstractMole{M,U} end

"""
    viscosity(T::Real,G::AbstractMole)

Laminar dynamic vicsosity `μ` is stress to normal acceleration ratio (Pa⋅s or slug⋅ft⁻¹⋅s⁻¹).
"""
@pure viscosity(::MoleGas{M,μ}) where {M,μ} = μ
@pure sutherlandviscosity(::MoleGas{M,μ,Tμ}) where {M,μ,Tμ} = Tμ

"""
    conductivity(T::Real,G::AbstractMole)

Laminar thermal conductivity `k` of temperature variation (W⋅m⁻¹⋅K⁻¹ or lb⋅s⁻¹⋅°R⁻¹).
"""
@pure conductivity(::MoleGas{M,μ,Tμ,k}) where {M,μ,Tμ,k} = k
@pure sutherlandconductivity(::MoleGas{M,μ,Tμ,k,Tk}) where {M,μ,Tμ,k,Tk} = Tk

gastext(::MoleGas) = "MoleGas{$(molarmass(G)),"

show(io::IO,G::MoleGas) = print(io,gastext(G),"μ=$(viscosity(G)),Tμ=$(sutherlandviscosity(G)),k=$(conductivity(G)),Tk=$(sutherlandconductivity(G))}")

for op ∈ (:viscosity,:conductivity)
    @eval @pure function $op(T::Real,G::MoleGas)
        Tk = $(Symbol(:sutherland,op))(G)
        ((2*$op(G))/sqrt(Tk))*(sqrt(T)/(1+Tk/T))
    end
end

@pure function viscond(μ0,Tμ,k0,Tk,T0=288.16)
    T = 2(T0^1.5)
    μ = μ0*sqrt(Tμ)*(T0+Tμ)/T
    k = k0*sqrt(Tk)*(T0+Tk)/T
    return μ,k
end
@pure function viscond(G,T=518.69)
    μ = viscosity(kelvin*T,G)*0.020886
    Tμ = rankine*sutherlandviscosity(G)
    k = conductivity(kelvin*T,G)*0.5778*778/3600
    Tk = rankine*sutherlandconductivity(G)
    return μ,Tμ,k,Tk
end

@pure heatratio(F::AbstractMole) = heatratio(288.16,F)
@pure heatvolume(F::AbstractMole) = heatvolume(288.16,F)
@pure heatpressure(F::AbstractMole) = heatpressure(288.16,F)

"""
    heatpressure(T::Real,G::AbstractMole) = heatvolume(T,G)+gasconstant(G)

Specific heat `cₚ` of ideal gas at constant pressure (J⋅kg⁻¹⋅K⁻¹ or ft⋅lb⋅slug⁻¹⋅°R⁻¹).
"""
@pure heatpressure(T::Real,G::MoleGas) = heatvolume(T,G)+gasconstant(G)

"""
    heatratio(T::Real,G::AbstractMole) = heatpressure(T,G)/heatvolume(T,G)

Specific heat ratio `γ` at constant pressure to constant volume of ideal gas (dimensionless).
"""
@pure heatratio(T::Real,G::MoleGas) = (gasconstant(G)/heatvolume(T,G))+1

"""
    energy(T::Real,G::AbstractMole) = heatvolume(T,G)*T

Specific energy `e` of ideal gas `enthalpy(T,G)-gasconstant(G)*T` (J⋅kg⁻¹ or ft⋅lb⋅slug⁻¹).
"""
@pure energy(T::Real,G::AbstractMole) = heatvolume(T,G)*T

"""
    enthalpy(T::Real,G::AbstractMole) = heatpressure(T,G)*T

Specific enthalpy `h` of ideal gas `energy(T,G)+gasconstant(G)*T` (J⋅kg⁻¹ or ft⋅lb⋅slug⁻¹).
"""
@pure enthalpy(T::Real,G::AbstractMole) = heatpressure(T,G)*T

"""
    freedom(T::Real,G::AbstractMole) = 2heatvolume(T,G)/gasconstant(G)

Degrees of freedom `f` with storage energy per mole `gasconstant(G)*T/2` (dimensionless).
"""
@pure freedom(T::Real,G::AbstractMole) = heatvolume(T,G)*(2/gasconstant(G))

"""
    prandtl(T::Real,G::AbstractMole) = viscosity(T,G)*heatpressure(T,G)/conductivity(T,G)

Prandtl number is the ratio of momentum diffusivity to heat diffusivity (dimensionless).
"""
@pure prandtl(T::Real,G::AbstractMole) = viscosity(T,G)*heatpressure(T,G)/conductivity(T,G)

"""
    sonicspeed(T::Real,G::AbstractMole) = sqrt(gasconstant(G)*heatratio(T,G)*T)

Speed of sound wave disturbance at isentropic conditions in ideal gas (m⋅s⁻¹ or ft⋅s⁻¹).
"""
@pure sonicspeed(T::Real,G::AbstractMole) = sqrt((gasconstant(G)*heatratio(T,G))*T)

export AtomicGas, DiatomicGas, TriatomicGas, SutherlandGas
export wavenumber, frequency, vibration

struct AtomicGas{M,μ,Tμ,k,Tk,U} <: MoleGas{M,μ,Tμ,k,Tk,U} end
struct DiatomicGas{M,ν,μ,Tμ,k,Tk,U} <: MoleGas{M,μ,Tμ,k,Tk,U} end
struct TriatomicGas{M,ν1,ν2,μ,Tμ,k,Tk,U} <: MoleGas{M,μ,Tμ,k,Tk,U} end
struct PentatomicGas{M,μ,Tμ,k,Tk,U} <: MoleGas{M,μ,Tμ,k,Tk,U} end
struct SutherlandGas{M,cᵥ,μ,Tμ,k,Tk,U} <: MoleGas{M,μ,Tμ,k,Tk,U} end

function AtomicGas(M,μ0,Tμ,k0,Tk,T0=288.16,U=Metric)
    μ,k = viscond(μ0,Tμ,k0,Tk,T0)
    AtomicGas{M,μ,Tμ,k,Tk,U}()
end
function DiatomicGas(M,ν,μ0,Tμ,k0,Tk,T0=288.16,U=Metric)
    μ,k = viscond(μ0,Tμ,k0,Tk,T0)
    DiatomicGas{M,ν,μ,Tμ,k,Tk,U}()
end
function TriatomicGas(M,ν1,ν2,μ0,Tμ,k0,Tk,T0=288.16,U=Metric)
    μ,k = viscond(μ0,Tμ,k0,Tk,T0)
    TriatomicGas{M,ν1,ν2,μ,Tμ,k,Tk,U}()
end
function PentatomicGas(M,μ0,Tμ,k0,Tk,T0=288.16,U=Metric)
    μ,k = viscond(μ0,Tμ,k0,Tk,T0)
    PentatomicGas{M,μ,Tμ,k,Tk,U}()
end
function SutherlandGas(M,cᵥ,μ0,Tμ,k0,Tk,T0,U=Metric)
    μ,k = viscond(μ0,Tμ,k0,Tk,T0)
    SutherlandGas{M,cᵥ,μ,Tμ,k,Tk,U}()
end

function AtomicGas(G::AtomicGas,T=518.69)
    μ,Tμ,k,Tk = viscond(G,T)
    AtomicGas(molarmass(G),μ,Tμ,k,Tk,T,English)
end
function DiatomicGas(G::DiatomicGas,T=518.69)
    μ,Tμ,k,Tk = viscond(G,T)
    M,ν = molarmass(G),meters(wavenumber(G))
    DiatomicGas(M,ν,μ,Tμ,k,Tk,T,English)
end
function TriatomicGas(G::TriatomicGas,T=518.69)
    μ,Tμ,k,Tk = viscond(G,T)
    ν1,ν2 = meters.(wavenumber(G))
    TriatomicGas(molarmass(G),ν1,ν2,μ,Tμ,k,Tk,T,English)
end
function PentatomicGas(G::PentatomicGas,T=518.69)
    μ,Tμ,k,Tk = viscond(G,T)
    PentatomicGas(molarmass(G),μ,Tμ,k,Tk,T,English)
end
function SutherlandGas(G::SutherlandGas,cᵥ,T=518.69)
    μ,Tμ,k,Tk = viscond(G,T)
    SutherlandGas(molarmass(G),cᵥ,μ,Tμ,k,Tk,T,English)
end

"""
    wavenumber(G::MoleGas) = frequency(G)/lightspeed(G)

Spectroscopic vibrational wavenumbers of polyatomic molecules if applicable (m⁻¹ or ft⁻¹).
"""
@pure wavenumber(::DiatomicGas{M,ν}) where {M,ν} = ν
@pure wavenumber(::TriatomicGas{M,ν1,ν2}) where {M,ν1,ν2} = ν1,ν2

"""
    wavelength(G::MoleGas) = 1/wavenumber(G)

Spectroscopic vibrational wavelength `λ` of polyatomic molecules if applicable (m or ft).
"""
@pure wavelength(G::MoleGas) = inv.(wavenumber(G))

"""
    frequency(G::MoleGas) = wavenumber(G)*lightspeed(G)

Spectroscopic vibrational frequencies `ν` of polyatomic molecules if applicable (Hz or s⁻¹).
"""
@pure frequency(G::MoleGas) = wavenumber(G).*lightspeed(G)

"""
    vibration(G::MoleGas) = frequency(G)*planck(G)/boltzmann(G)

Vibrational temperature `θ` of polyatomic molecules if applicable (K or °R).
"""
@pure vibration(G::MoleGas) = frequency(G).*(planck(G)/boltzmann(G)/1.2)
@pure vibration(θT::Float64) = (eθT = exp(θT); θT^2*eθT/(eθT-1)^2)

"""
    heatvolume(T::Real,G::AbstractMole) = translation + rotation + vibration

Specific heat `cᵥ` of ideal gas at constant volume (J⋅kg⁻¹⋅K⁻¹ or ft⋅lb⋅slug⁻¹⋅°R⁻¹).
"""
@pure heatvolume(::Real,G::AtomicGas) = (3/2)*gasconstant(G)
@pure heatvolume(::Real,::SutherlandGas{M,cᵥ}) where {M,cᵥ} = cᵥ
@pure function heatvolume(T::Real,G::DiatomicGas)
    R,θ = gasconstant(G),vibration(G)
    R*((5/2)+vibration(θ/T))
end
@pure function heatvolume(T::Real,G::TriatomicGas)
    R = gasconstant(G)
    θ1,θ2 = vibration(G)
    R*((5/2)+vibration(θ1/T)+vibration(θ2/T))
end

gastext(G::AtomicGas) = "AtomicGas{M=$(molarmass(G)),"
gastext(G::DiatomicGas) = "DiatomicGas{M=$(molarmass(G)),ν=$(vibration(G)),"
gastext(G::TriatomicGas) = "TriatomicGas{M=$(molarmass(G)),ν₁=$(vibration(G)[1]),ν₂=$(vibration(G)[2]),"
gastext(G::SutherlandGas) = "Gas{M=$(molarmass(G)),cᵥ=$(heatvolume(G)),cₚ=$(heatpressure(G)),"

struct Mixture{M,N,C,U} <: AbstractMole{M,U}
    f::Values{N,Float64}
end

Mixture{N,C,U}(f) where {N,C,U} = Mixture{f⋅molarmass.(C),N,C,U}(f)
Mixture{C,U}(f) where {C,U} = Mixture{length(C),C,U}(f)
Mixture{C}(f) where C = Mixture{C,Metric}(f)

@pure chemical(::Mixture{M,N,C}) where {M,N,C} = C
@pure molecules(::Mixture{M,N,C}) where {M,N,C} = Values(C)
@pure fractions(M::Mixture) = M.f
@pure Base.length(::Mixture{M,N}) where {M,N} = N
@pure Base.getindex(M::Mixture,i::Int) = chemical(M)[i]
@pure heatratio(T::Real,M::Mixture) = heatpressure(T,M)/heatvolume(T,M)

for op ∈ (:viscosity,:conductivity,:heatvolume,:heatpressure)
    @eval $op(T::Real,M::Mixture) = fractions(M)⋅Values($op.(T,chemical(M)))
end
for op ∈ (:viscosity,:conductivity,:sutherlandviscosity,:sutherlandconductivity)
    @eval $op(M::Mixture) = fractions(M)⋅Values($op.(chemical(M)))
end

Base.:*(f::Float64,G::AbstractMole{M,U}) where {M,U} = Mixture{M,1,(G,),U}(Values(f))

function Base.:+(m::Mixture...)
    N = sum(length.(m))
    Mixture{N,Tuple(vcat(molecules.(m)...)),units(m[1])}(Values{N}(vcat(fractions.(m)...)))
end

