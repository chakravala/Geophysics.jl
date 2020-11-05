module Geophysics

#   This file is part of Geophysics.jl. It is licensed under the MIT license
#   Geophysics Copyright (C) 2020 Michael Reed

using AbstractTensors
import Base: @pure, show, display

export Fluid, Gas, FluidState, Atmosphere, Weather
export heatvolume, heatpressure, heatratio, molarmass, gasconstant
export fluid, radius, gravity, layer, sutherlandviscosity, sutherlandconductivity

# chemistry

const Boltzmann,Avogadro = 1.380649e-23,6.02214076e23
const Rᵤ,RᵤEnglish = Boltzmann*Avogadro,49.72
const AvogadroEnglish = 453.59237Avogadro
const BoltzmannEnglish = RᵤEnglish/AvogadroEnglish
const MEnglish,kEnglish = Rᵤ/RᵤEnglish,0.5778

# fluid models

abstract type Fluid{M,cᵥ,cₚ} end

@pure heatvolume(::Fluid{M,cᵥ}) where {M,cᵥ} = cᵥ
@pure heatpressure(::Fluid{M,cᵥ,cₚ}) where {M,cᵥ,cₚ} = cₚ
@pure heatratio(F::Fluid) = heatpressure(F)/heatvolume(F)
@pure molarmass(::Fluid{M}) where M = M
@pure gasconstant(F::Fluid) = Rᵤ/molarmass(F)

struct Gas{M,cᵥ,cₚ,μ,Tμ,k,Tk} <: Fluid{M,cᵥ,cₚ} end
struct Liquid{M,cᵥ,cₚ} <: Fluid{M,cᵥ,cₚ} end

show(io::IO,F::Fluid) = print(io,"Fluid{M=$(molarmass(F)),cᵥ=$(heatvolume(F)),cₚ=$(heatpressure(F))}")
show(io::IO,F::Gas) = print(io,"Gas{M=$(molarmass(F)),cᵥ=$(heatvolume(F)),cₚ=$(heatpressure(F)),μ=$(viscosity(F)),Tμ=$(sutherlandviscosity(F)),k=$(conductivity(F)),Tk=$(sutherlandconductivity(F))}")

function Gas(M,cᵥ,cₚ,μ0,Tμ,k0,Tk,T0)
    T = 2(T0^1.5)
    μ = μ0*sqrt(Tμ)*(T0+Tμ)/T
    k = k0*sqrt(Tk)*(T0+Tk)/T
    Gas{M,cᵥ,cₚ,μ,Tμ,k,Tk}()
end

@pure function conductivity(T::Real,G::Gas)
    Tk = sutherlandconductivity(G)
    (2conductivity(G)/sqrt(Tk))*(sqrt(T)/(1+Tk/T))
end
@pure function viscosity(T::Real,G::Gas)
    Tμ = sutherlandviscosity(G)
    (2viscosity(G)/sqrt(Tμ))*(sqrt(T)/(1+Tμ/T))
end

@pure viscosity(F::Gas{M,cᵥ,cₚ,μ}) where {M,cᵥ,cₚ,μ} = μ
@pure conductivity(F::Gas{M,cᵥ,cₚ,μ,Tμ,k}) where {M,cᵥ,cₚ,μ,Tμ,k} = k
@pure sutherlandviscosity(::Gas{M,cᵥ,cₚ,μ,Tμ}) where {M,cᵥ,cₚ,μ,Tμ} = Tμ
@pure sutherlandconductivity(::Gas{M,cᵥ,cₚ,μ,Tμ,k,Tk}) where {M,cᵥ,cₚ,μ,Tμ,k,Tk} = Tk

struct FluidState{f}
    T::Float64
    P::Float64
end

(G::Gas)(T=288.15,P=101325.0) = FluidState{G}(T,P)
(L::Liquid)(T,P) = FluidState{L}(T,P)

@pure fluid(::FluidState{f}) where f = f
@pure heatpressure(F::FluidState) = heatpressure(fluid(F))
@pure heatvolume(F::FluidState) = heatvolume(fluid(F))
@pure heatratio(F::FluidState) = heatratio(fluid(F))
@pure molarmass(F::FluidState) = molarmass(fluid(F))
@pure gasconstant(F::FluidState) = gasconstant(fluid(F))

@pure temperature(F::FluidState) = F.T
@pure pressure(F::FluidState) = F.P
@pure density(F::FluidState) = (pressure(F)/temperature(F))/gasconstant(F)
@pure conductivity(F::FluidState) = conductivity(temperature(F),fluid(F))
@pure viscosity(F::FluidState) = viscosity(temperature(F),fluid(F))
@pure kinematic(F::FluidState) = viscosity(F)/density(F)
@pure volume(F::FluidState) = inv(density(F))
@pure energy(F::FluidState) = heatvolume(F)*temperature(F)
@pure enthalpy(F::FluidState) = heatpressure(F)*temperature(F)
@pure heatcapacity(F::FluidState) = heatpressure(F)*density(F)
@pure diffusivity(F::FluidState) = conductivity(F)/heatcapacity(F)
@pure elasticity(F::FluidState) = heatratio(F)*pressure(F)
@pure prandtl(F::FluidState) = heatpressure(F)*(viscosity(F)/conductivity(F))
@pure sonicspeed(F::FluidState) = sqrt(heatratio(F)*gasconstant(F))*sqrt(temperature(F))
@pure impedance(F::FluidState) = density(F)*sonicspeed(F)

# global geophysics

struct Atmosphere{f,n}
    a::Values{n,Float64} # lapse rate
    h::Values{n,Float64} # altitude
    Atmosphere{f}(a::Values{n},h::Values{n}) where {f,n} = new{f,n}(a,h)
end

function display(A::Atmosphere)
    println(typeof(A))
    println(" a = ",A.a)
    println(" h = ",A.h)
end

# local weather models

struct Weather{r,g,f,n} # radius, acceleration
    A::Atmosphere{f,n} # altitude lapse rate
    T::Values{n,Float64} # temperature
    p::Values{n,Float64} # pressure
    ρ::Values{n,Float64} # density
    Weather{r,g}(A::Atmosphere{f},T::Values{n},p::Values{n},ρ::Values{n}) where {r,g,f,n} = new{r,g,f,n}(A,T,p,ρ)
end

function Weather{r,g}(A::Atmosphere{f,n},F::FluidState) where {r,g,f,n}
    T = zeros(Variables{n,Float64})
    p = zeros(Variables{n,Float64})
    ρ = zeros(Variables{n,Float64})
    T0,p0 = temperature(F),pressure(F)
    R,γ = gasconstant(F),heatratio(F)
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
    Weather{r,g}(A,Values(T),Values(p),Values(ρ))
end

(A::Atmosphere)(F::FluidState=Air(288.16),r=6.356766e6,g=9.80665) = Weather{r,g}(A,F)
(A::Atmosphere)(T,p=101325.0,r=6.356766e6,g=9.80665) = Weather{r,g}(A,fluid(A)(T,p))
(W::Weather)(hG::Real=0) = W(hG,layer(hG,W))
function (W::Weather)(hG::Real,i)
    g,R = gravity(W),gasconstant(W)
    T0,a,h0,p = W[i]; h = altgeopotent(hG,W); T = _temperature(h,i,W)
    FluidState{fluid(W)}(T,if iszero(a)
        p*exp((-g/R)*(h-h0)/T)
    else
        p*(T/T0)^((-g/R)/a)
    end)
end

function display(W::Weather)
    println(typeof(W))
    println(" a = ",W.A.a)
    println(" h = ",W.A.h)
    println(" T = ",W.T)
    println(" P = ",W.p)
    println(" ρ = ",W.ρ)
end

@pure @inline Base.getindex(W::Weather,i::Int) = W.T[i],W.A.a[i],W.A.h[i],W.p[i],W.ρ[i]
@pure @inline Base.getindex(W::Weather,::Val{i}) where i = getindex(W,i)
@pure @inline lapserate(h::Real,W::Weather=Standard) = W[layer(h,W)]
@pure layer(h::Real,W::Weather=Standard) = h≤W.A.h[1] ? 1 : (i=findfirst(x->x≥h,W.A.h); isnothing(i) ? length(W.A.h) : i-1)

@pure fluid(::Atmosphere{f}) where f = f
@pure fluid(::Weather{r,g,f}=Standard) where {r,g,f} = f
@pure heatvolume(W::Weather=Standard) = heatvolume(fluid(W))
@pure heatpressure(W::Weather=Standard) = heatpressure(fluid(W))
@pure heatratio(W::Weather=Standard) = heatratio(fluid(W))
@pure molarmass(W::Weather=Standard) = molarmass(fluid(W))
@pure gasconstant(W::Weather=Standard) = gasconstant(fluid(W))
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

@pure conductivity(hG::Real,i,W::Weather=Standard) = conductivity(temperature(hG,i,W),fluid(W))

# μ = viscosity

@pure viscosity(hG::Real,i,W::Weather=Standard) = viscosity(temperature(hG,i,W),fluid(W))

# ν = kinematic viscosity

@pure function kinematic(hG::Real,i,W::Weather=Standard)
    g,R = gravity(W),gasconstant(W)
    T0,a,h0,_,ρ = W[i]; h = altgeopotent(hG,W); T = _temperature(h,i,W)
    viscosity(T,fluid(W))/(if iszero(a)
        ρ*exp((-g/R)*(h-h0)/T)
    else
        ρ*(T/T0)^((-g/R)/a-1)
    end)
end

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

@pure function diffusivity(hG::Real,i,W::Weather=Standard)
    g,R = gravity(W),gasconstant(W)
    T0,a,h0,_,ρ = W[i]; h = altgeopotent(hG,W); T = _temperature(h,i,W)
    conductivity(T,fluid(W))/(if iszero(a)
        ρ*exp((-g/R)*(h-h0)/T)
    else
        ρ*(T/T0)^((-g/R)/a-1)
    end)/heatpressure(W)
end

# bulk elasticity

@pure elasticity(hG,i,W::Weather=Standard) = heatratio(W)*pressure(hG,i,W)

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

@pure function impedance(hG::Real,i,W::Weather=Standard)
    g,R = gravity(W),gasconstant(W)
    T0,a,h0,_,ρ = W[i]; h = altgeopotent(hG,W); T = _temperature(h,i,W)
    ρ*(if iszero(a)
        exp((-g/R)*(h-h0)/T)
    else
        (T/T0)^((-g/R)/a-1)
    end)*sqrt(heatratio(W)*R)*sqrt(T)
end

# common interface

for op ∈ (:temperature,:pressure,:density,:sonicspeed,:weight,:volume,:potential,:impedance,:grashof,:prandtl,:diffusivity,:heatcapacity,:enthalpy,:energy,:kinematic,:viscosity,:conductivity,:elasticity)
    opratio = Symbol(op,:ratio)
    @eval begin
        export $op
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
