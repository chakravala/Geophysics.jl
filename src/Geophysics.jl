module Geophysics

#   This file is part of Geophysics.jl. It is licensed under the AGPL license
#   Geophysics Copyright (C) 2020 Michael Reed


using AbstractTensors, LinearAlgebra
import Base: @pure, show, display

export FluidState, Atmosphere, Weather

export fluid, radius, gravity, layer

# thermodynamics

include("chemistry.jl")

struct FluidState{f}
    T::Float64
    P::Float64
end

for Gas ∈ (:SutherlandGas,:Mixture,:AtomicGas,:DiatomicGas,:TriatomicGas)
    @eval (G::$Gas)(T=288.15,P=101325.0) = FluidState{G}(T,P)
end

@pure fluid(::FluidState{f}) where f = f
@pure temperature(F::FluidState) = F.T
@pure pressure(F::FluidState) = F.P

for op ∈ Properties
    @eval @pure $op(F::FluidState) = $op(fluid(F))
end
for op ∈ Intrinsic
    @eval @pure $op(F::FluidState) = $op(temperature(F),fluid(F))
end

@pure density(F::FluidState) = (pressure(F)/temperature(F))/gasconstant(F)
@pure kinematic(F::FluidState) = viscosity(F)/density(F)
@pure volume(F::FluidState) = inv(density(F))
@pure heatcapacity(F::FluidState) = heatpressure(F)*density(F)
@pure diffusivity(F::FluidState) = conductivity(F)/heatcapacity(F)
@pure elasticity(F::FluidState) = heatratio(F)*pressure(F)
@pure impedance(F::FluidState) = density(F)*sonicspeed(F)

# global geophysics

struct Atmosphere{f,n}
    a::Values{n,Float64} # lapse rate
    h::Values{n,Float64} # altitude
    m::Values{n,Float64} # molar rate
    Atmosphere{f}(a::Values{n},h::Values{n},m::Values{n}) where {f,n} = new{f,n}(a,h,m)
end

Atmosphere{f}(a::Values{n},h::Values{n}) where {f,n} = Atmosphere{f}(a,h,zeros(Values{n}))

function display(A::Atmosphere)
    println(typeof(A))
    println(" a = ",A.a)
    println(" h = ",A.h)
    println(" m = ",A.m)
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
@pure radius(::Weather{r}=Standard) where r = r
@pure gravity(::Weather{r,g}=Standard) where {r,g} = g

for op ∈ Properties
    @eval @pure $op(W::Weather=Standard) = $op(fluid(W))
end

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

# k = thermal conductivity
# μ = viscosity
# cᵥ = specific heat (constant volume)
# cₚ = specific heat (constant pressure)
# γ = specific heat ratio
# Pr = Prandtl number
# a = speed of sound
# e = energy
# h = enthalpy

for op ∈ Intrinsic
    @eval @pure $op(hG::Real,i,W::Weather=Standard) = $op(temperature(hG,i,W),fluid(W))
end

# p = pressure

@pure function pressure(h::Real,T,i,W::Weather=Standard)
    g,R = gravity(W),gasconstant(W)
    T0,a,h0,p = W[i]
    p*(if iszero(a)
        exp((-g/R)*(h-h0)/T)
    else
        (T/T0)^((-g/R)/a)
    end)
end
@pure function pressure(hG::Real,i,W::Weather=Standard)
    h = altgeopotent(hG,W)
    pressure(h,_temperature(h,i,W),i,W)
end

# ρ = density

@pure function density(h::Real,T,i,W::Weather=Standard)
    g,R = gravity(W),gasconstant(W)
    T0,a,h0,_,ρ = W[i]
    ρ*(if iszero(a)
        exp((-g/R)*(h-h0)/T)
    else
        (T/T0)^((-g/R)/a-1)
    end)
end
@pure function density(hG::Real,i,W::Weather=Standard)
    h = altgeopotent(hG,W)
    density(h,_temperature(h,i,W),i,W)
end

# ν = kinematic viscosity

@pure kinematic(h::Real,T,i,W::Weather=Standard) = viscosity(T,fluid(W))/density(h,T,i,W)
@pure function kinematic(hG::Real,i,W::Weather=Standard)
    h = altgeopotent(hG,W)
    kinematic(h,_temperature(h,i,W),i,W)
end

# heat capacity

@pure function heatcapacity(hG::Real,i,W::Weather=Standard)
    h = altgeopotent(hG,W)
    T = _temperature(h,i,W)
    heatpressure(T,fluid(W))*density(h,T,i,W)
end

# thermal diffusivity

@pure function diffusivity(hG::Real,i,W::Weather=Standard)
    F = fluid(W)
    h = altgeopotent(hG,W)
    T = _temperature(h,i,W)
    conductivity(T,F)/heatpressure(T,F)/density(h,T,i,W)
end

# bulk elasticity

@pure function elasticity(hG,i,W::Weather=Standard)
    h = altgeopotent(hG,W)
    T = _temperature(h,i,W)
    heatratio(T,fluid(W))*pressure(h,T,i,W)
end

# characteristic specific acoustic impedance

@pure function impedance(hG::Real,i,W::Weather=Standard)
    h = altgeopotent(hG,W)
    T = _temperature(h,i,W)
    density(h,T,i,W)*sonicspeed(T,fluid(W))
end

# Grashof number

@pure function grashof(hG::Real,i,W::Weather=Standard)
    h = altgeopotent(hG,W)
    T = _temperature(h,i,W)
    gravity(hG,W)*(temperature(W)-T)*(hG^3)/(T*kinematic(h,T,i,W)^2)
end

# specific weight

@pure weight(hG::Real,i,W::Weather=Standard) = density(hG,i,W)*gravity(hG,W)

# specific volume

@pure volume(hG::Real,i,W::Weather=Standard) = inv(density(hG,i,W))

# specific potential energy

@pure potential(hG::Real,i,W::Weather=Standard) = gravity(hG,W)*hG

# common interface

for op ∈ (:temperature,:pressure,:density,:weight,:volume,:potential,:impedance,:grashof,:diffusivity,:heatcapacity,:kinematic,:elasticity,Intrinsic...)
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
