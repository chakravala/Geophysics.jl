module Geophysics

#   This file is part of Geophysics.jl. It is licensed under the MIT license
#   Geophysics Copyright (C) 2020 Michael Reed

using AbstractTensors
import Base: @pure

export Fluid, Gas, Atmosphere, Weather
export heatvolume, heatpressure, heatratio, molarmass, gasconstant
export fluid, radius, gravity, layer, sutherlandviscosity, sutherlandconductivity

# fluid models

abstract type Fluid{M,cᵥ,cₚ} end
struct Gas{M,cᵥ,cₚ,Tμ,Tk} <: Fluid{M,cᵥ,cₚ} end
struct Liquid{M,cᵥ,cₚ} <: Fluid{M,cᵥ,cₚ} end

const Boltzmann,Avogadro = 1.380649e-23,6.02214076e23
const Rᵤ,RᵤEnglish = Boltzmann*Avogadro,49.72
const AvogadroEnglish = 453.59237Avogadro
const BoltzmannEnglish = RᵤEnglish/AvogadroEnglish
const MEnglish,kEnglish = Rᵤ/RᵤEnglish,0.5778

@pure heatvolume(::Fluid{M,cᵥ}) where {M,cᵥ} = cᵥ
@pure heatpressure(::Fluid{M,cᵥ,cₚ}) where {M,cᵥ,cₚ} = cₚ
@pure heatratio(F::Fluid) = heatpressure(F)/heatvolume(F)
@pure molarmass(::Fluid{M}) where M = M
@pure gasconstant(F::Fluid) = Rᵤ/molarmass(F)
@pure sutherlandviscosity(::Gas{M,cᵥ,cₚ,Tμ}) where {M,cᵥ,cₚ,Tμ} = Tμ
@pure sutherlandconductivity(::Gas{M,cᵥ,cₚ,Tμ,Tk}) where {M,cᵥ,cₚ,Tμ,Tk} = Tk

struct Atmosphere{f,n}
    a::Values{n,Float64} # lapse rate
    h::Values{n,Float64} # altitude
    Atmosphere{f}(a::Values{n},h::Values{n}) where {f,n} = new{f,n}(a,h)
end

# local weather models

struct Weather{r,g,f,n,μ,k} # radius, acceleration
    A::Atmosphere{f,n} # altitude lapse rate
    T::Values{n,Float64} # temperature
    p::Values{n,Float64} # pressure
    ρ::Values{n,Float64} # density
    Weather{r,g,μ,k}(A::Atmosphere{f},T::Values{n},p::Values{n},ρ::Values{n}) where {r,g,μ,k,f,n} = new{r,g,f,n,μ,k}(A,T,p,ρ)
end

function Weather{r,g,μ,k}(T0,p0,A::Atmosphere{f,n}) where {r,g,μ,k,f,n} T = zeros(Variables{n,Float64})
    p = zeros(Variables{n,Float64})
    ρ = zeros(Variables{n,Float64})
    R,γ = gasconstant(f),heatratio(f)
    T[1] = T0; p[1] = p0; ρ[1] = p0/(R*T0)
    for i ∈ 2:n
        Δh = A.h[i]-A.h[i-1]
        T[i] = T[i-1]+A.a[i-1]*Δh
        vp,vρ = if iszero(A.a[i-1])
            val = exp((-g/R)*Δh/T[i])
            val,val
        else
            t,gRa = T[i]/T[i-1],(-g/R)/A.a[i-1]
            t^gRa,t^(gRa-1)
        end
        p[i] = p[i-1]*vp
        ρ[i] = ρ[i-1]*vρ
    end
    Weather{r,g,μ,k}(A,Values(T),Values(p),Values(ρ))
end

(A::Atmosphere)(r,g,μ,k,T,p) = Weather{r,g,μ,k}(T,p,A)

@pure @inline Base.getindex(W::Weather,i::Int) = W.T[i],W.A.a[i],W.A.h[i],W.p[i],W.ρ[i]
@pure @inline Base.getindex(W::Weather,::Val{i}) where i = getindex(W,i)
@pure @inline lapserate(h::Real,W::Weather=Standard) = W[layer(h,W)]
@pure layer(h::Real,W::Weather=Standard) = h≤W.A.h[1] ? 1 : (i=findfirst(x->x≥h,W.A.h); isnothing(i) ? length(W.A.h) : i-1)

@pure fluid(::Weather{r,g,f}=Standard) where {r,g,f} = f
@pure heatvolume(W::Weather=Standard) = heatvolume(fluid(W))
@pure heatpressure(W::Weather=Standard) = heatpressure(fluid(W))
@pure heatratio(W::Weather=Standard) = heatratio(fluid(W))
@pure gasconstant(W::Weather=Standard) = gasconstant(fluid(W))
@pure viscosity(::Weather{r,g,f,n,μ}=Standard) where {r,g,f,n,μ} = μ
@pure conductivity(::Weather{r,g,f,n,μ,k}=Standard) where {r,g,f,n,μ,k} = k
@pure radius(::Weather{r}=Standard) where r = r
@pure gravity(::Weather{r,g}=Standard) where {r,g} = g

# h = geopotential altitude, hG = geometric altitude

@pure @inline altabs(hG::Real=0,W::Weather=Standard) = radius(W)+hG
@pure @inline altratio(hG::Real,W::Weather=Standard) = radius(W)/altabs(hG,W)
@pure altgeopotent(hG::Real,W::Weather=Standard) = hG*altratio(hG,W)
@pure altgeometric(h::Real,W::Weather=Standard) = (r=radius(W); r/(r/h-1))
@pure gravity(hG::Real,W::Weather=Standard) = (gravity(W)*radius(W)^2)/altabs(hG,W)^2

# T = temperature

@pure temperature(h::Real,i,W::Weather=Standard) = _temperature(altgeopotent(h,W),i,W)
@pure function _temperature(h::Real,i,W::Weather=Standard)
    T0,a0,h0 = W[i]
    return T0+a0*(h-h0)
end

# p = pressure

@pure function pressure(hG::Real,i,W::Weather=Standard)
    g,R = gravity(W),gasconstant(W)
    T0,a,h0,p = W[i]; h = altgeopotent(hG,W); T = _temperature(h,i,W)
    if iszero(a)
        p*exp((-g/R)*(h-h0)/T)
    else
        p*(T/T0)^((-g/R)/a)
    end
end

# ρ = density

@pure function density(hG::Real,i,W::Weather=Standard)
    g,R = gravity(W),gasconstant(W)
    T0,a,h0,_,ρ = W[i]; h = altgeopotent(hG,W); T = _temperature(h,i,W)
    if iszero(a)
        ρ*exp((-g/R)*(h-h0)/T)
    else
        ρ*(T/T0)^((-g/R)/a-1)
    end
end

# k = thermal conductivity

@pure conductivity(hG::Real,W::Weather=Standard) = conductivity(hG,layer(hG,W),W)
@pure function conductivity(hG::Real,i,W::Weather=Standard,G::Gas=fluid(W))
    T,T0,S = temperature(hG,i,W),temperature(W),sutherlandconductivity(G)
    return (conductivity(W)*(T0+S)/T0^1.5)*(sqrt(T)/(1+S/T))
end

# μ = viscosity

@pure viscosity(hG::Real,W::Weather=Standard) = viscosity(hG,layer(hG,W),W,fluid(W))
@pure function viscosity(hG::Real,i,W::Weather=Standard,G::Gas=fluid(W))
    T,T0,S = temperature(hG,i,W),temperature(W),sutherlandviscosity(G)
    return (viscosity(W)*(T0+S)/T0^1.5)*(sqrt(T)/(1+S/T))
end

# ν = kinematic viscosity

@pure kinematic(hG::Real,i,W::Weather=Standard) = viscosity(hG,i,W)/density(hG,i,W)

# specific weight

@pure weight(hG::Real,i,W::Weather=Standard) = density(hG,i,W)*gravity(hG,W)

# specific volume

@pure volume(hG::Real,i,W::Weather=Standard) = inv(density(hG,i,W))

# specific potential energy

@pure potential(hG::Real,i,W::Weather=Standard) = gravity(hG,W)*hG

# e = energy

@pure energy(hG::Real,i,W::Weather=Standard) = heatvolume(W)*temperature(hG,i,W)

# h = enthalpy

@pure enthalpy(hG::Real,i,W::Weather=Standard) = heatpressure(W)*temperature(hG,i,W)

# heat capacity

@pure heatcapacity(hG::Real,i,W::Weather=Standard) = heatpressure(W)*density(hG,i,W)

# thermal diffusivity

@pure diffusivity(hG::Real,i,W::Weather=Standard) = conductivity(hG,i,W)/heatcapacity(hG,i,W)

# Prandtl number

@pure prandtl(hG::Real,i,W::Weather=Standard) = viscosity(hG,i,W)*heatpressure(W)/conductivity(hG,i,W)

# Grashof number

@pure function grashof(hG::Real,i,W::Weather=Standard)
    T = temperature(hG,i,W)
    gravity(hG,W)*(temperature(W)-T)*(hG^3)/(T*kinematic(hG,i,W)^2)
end

# speed of sound

@pure sonicspeed(hG::Real,i,W::Weather=Standard) = sqrt(heatratio(W)*gasconstant(W))*sqrt(temperature(hG,i,W))

# characteristic specific acoustic impedance

@pure impedance(hG::Real,i,W::Weather=Standard) = density(hG,i,W)*sonicspeed(hG,i,W)

# common interface

for op ∈ (:temperature,:pressure,:density,:sonicspeed,:weight,:volume,:potential,:impedance,:grashof,:prandtl,:diffusivity,:heatcapacity,:enthalpy,:energy,:kinematic,:viscosity,:conductivity)
    opratio = Symbol(op,:ratio)
    @eval export $op
    op ∉ (:viscosity,:conductivity) && @eval begin
        @pure $op(W::Weather=Standard) = $op(0,W)
        @pure $op(hG::Real,W::Weather=Standard) = $op(hG,layer(hG,W),W)
    end
    op ∉ (:heatcapacity,:grashof,:potential) && @eval begin
        export $opratio
        @pure $opratio(hG::Real,W::Weather=Standard) = $opratio(hG,layer(hG,W),W)
        @pure $opratio(hG::Real,i,W::Weather=Standard) = $op(hG,i,W)/$op(W)
    end
end

# data

include("planets.jl")

end # module
