
#   This file is part of Geophysics.jl. It is licensed under the MIT license
#   Geophysics Copyright (C) 2020 Michael Reed

export rankine, kelvin, Metric, English

# chemistry units

struct UnitSystem{k,N,h,c} end
@pure boltzmann(::UnitSystem{k}) where k = k
@pure avogadro(::UnitSystem{k,N}) where {k,N} = N
@pure planck(::UnitSystem{k,N,h}) where {k,N,h} = h
@pure lightspeed(::UnitSystem{k,N,h,c}) where {k,N,h,c} = c
@pure universal(U::UnitSystem) = boltzmann(U)*avogadro(U)

const Planck = 6.62607015e-34
const PlanckEnglish = Planck/1055.056
const Boltzmann,Avogadro = 1.380649e-23,6.02214076e23
const AvogadroEnglish = 453.59237Avogadro
const BoltzmannEnglish = 49.72/AvogadroEnglish
const Metric = UnitSystem{Boltzmann,Avogadro,Planck,299792458.0}()
const English = UnitSystem{BoltzmannEnglish,AvogadroEnglish,PlanckEnglish,983569553.8}()

@pure rankine(T) = (9/5)T
@pure kelvin(T) = (5/9)T

const Constants = (:boltzmann,:avogadro,:planck,:lightspeed,:universal)

const Properties = (:molarmass,:gasconstant,Constants...)

const Intrinsic = (:viscosity,:conductivity,:heatvolume,:heatpressure,:heatratio,:prandtl,:sonicspeed,:freedom,:energy,:enthalpy)

# fluid models

export molarmass, units, gasconstant
export sutherlandviscosity, sutherlandconductivity

abstract type AbstractMole{M,U} end
@pure molarmass(::AbstractMole{M}) where M = M
@pure units(::AbstractMole{M,U}) where {M,U} = U
@pure gasconstant(F::AbstractMole) = universal(F)/molarmass(F)

for op ∈ Constants
    @eval begin; export $op
        @pure $op(F::AbstractMole) = $op(units(F))
    end
end

abstract type MoleGas{M,μ,Tμ,k,Tk,U} <: AbstractMole{M,U} end
@pure viscosity(::MoleGas{M,μ}) where {M,μ} = μ
@pure conductivity(::MoleGas{M,μ,Tμ,k}) where {M,μ,Tμ,k} = k
@pure sutherlandviscosity(::MoleGas{M,μ,Tμ}) where {M,μ,Tμ} = Tμ
@pure sutherlandconductivity(::MoleGas{M,μ,Tμ,k,Tk}) where {M,μ,Tμ,k,Tk} = Tk

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
    μ = viscosity(kelvin(T),G)*0.020886
    Tμ = rankine(sutherlandviscosity(G))
    k = conductivity(kelvin(T),G)*0.5778
    Tk = rankine(sutherlandconductivity(G))
    return μ,Tμ,k,Tk
end

@pure heatratio(F::AbstractMole) = heatratio(288.16,F)
@pure heatvolume(F::AbstractMole) = heatvolume(288.16,F)
@pure heatpressure(F::AbstractMole) = heatpressure(288.16,F)
@pure heatpressure(T::Real,G::MoleGas) = heatvolume(T,G)+gasconstant(G)
@pure heatratio(T::Real,G::MoleGas) = (gasconstant(G)/heatvolume(T,G))+1

@pure energy(T::Real,G::AbstractMole) = heatvolume(T,G)*T
@pure enthalpy(T::Real,G::AbstractMole) = heatpressure(T,G)*T
@pure freedom(T::Real,G::AbstractMole) = heatvolume(T,G)*(2/gasconstant(G))
@pure prandtl(T::Real,G::AbstractMole) = viscosity(T,G)*heatpressure(T,G)/conductivity(T,G)
@pure sonicspeed(T::Real,G::AbstractMole) = sqrt((gasconstant(G)*heatratio(T,G))*T)

gastext(::MoleGas) = "MoleGas{$(molarmass(G)),"

show(io::IO,G::MoleGas) = print(io,gastext(G),"μ=$(viscosity(G)),Tμ=$(sutherlandviscosity(G)),k=$(conductivity(G)),Tk=$(sutherlandconductivity(G))}")

export AtomicGas, DiatomicGas, TriatomicGas, SutherlandGas

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
    DiatomicGas{M,ν*lightspeed(U)/1.2,μ,Tμ,k,Tk,U}()
end
function TriatomicGas(M,ν1,ν2,μ0,Tμ,k0,Tk,T0=288.16,U=Metric)
    μ,k = viscond(μ0,Tμ,k0,Tk,T0)
    c = lightspeed(U)
    TriatomicGas{M,ν1*c/1.2,ν2*c/1.2,μ,Tμ,k,Tk,U}()
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
    M,ν = molarmass(G),30*vibration(G)/lightspeed(English)
    DiatomicGas(M,ν,μ,Tμ,k,Tk,T,English)
end
function TriatomicGas(G::TriatomicGas,T=518.69)
    μ,Tμ,k,Tk = viscond(G,T)
    ν1,ν2 = vibration(G)./(lightspeed(English)/30)
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

@pure vibration(θ::Float64) = (eθ = exp(θ); θ^2*eθ/(eθ-1)^2)
@pure vibration(::DiatomicGas{M,ν}) where {M,ν} = ν
@pure vibration(::TriatomicGas{M,ν1,ν2}) where {M,ν1,ν2} = ν1,ν2

@pure heatvolume(::Real,G::AtomicGas) = (3/2)*gasconstant(G)
@pure heatvolume(::Real,::SutherlandGas{M,cᵥ}) where {M,cᵥ} = cᵥ
@pure function heatvolume(T::Real,G::DiatomicGas)
    R,k,h,ν = gasconstant(G),boltzmann(G),planck(G),vibration(G)
    R*((5/2)+vibration((ν*(h/k))/T))
end
@pure function heatvolume(T::Real,G::TriatomicGas)
    R,k,h = gasconstant(G),boltzmann(G),planck(G)
    hk = h/k
    ν1,ν2 = vibration(G)
    R*((5/2)+vibration((ν1*hk)/T)+vibration((ν2*hk)/T))
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

