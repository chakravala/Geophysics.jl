module Geophysics

#   This file is part of Geophysics.jl. It is licensed under the AGPL license
#   Geophysics Copyright (C) 2020 Michael Reed


using AbstractTensors, LinearAlgebra
import Base: @pure, show, display

export FluidState, Atmosphere, Weather

export fluid, radius, gravity, layer

decibel(I₁,I₂=intensity()) = 10log10(I₁/I₂)
gage(P::Real,P0::Real=pressure()) = P-P0

# thermodynamics

include("units.jl")
include("chemistry.jl")

"""
    FluidState{f}

Thermodynamic state of fluid `f` at temperature `T` and pressure `P`.
Induces derived values `fluid`, `temperature`, `pressure`, `density`, `volume`, `kinematic`, `heatcapacity`, `diffusivity`, `elasticity`, `impedance`, `intensity`, and values associated with `f::AbstractMole` derivations.
"""
struct FluidState{f}
    T::Float64
    P::Float64
end

for Gas ∈ (:SutherlandGas,:Mixture,:AtomicGas,:DiatomicGas,:TriatomicGas)
    @eval (G::$Gas)(T=288.15,P=101325.0) = FluidState{G}(T,P)
end

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

"""
    pressure(F::FluidState)

Absolute force per unit area `P` exerted by `F` (Pa or slug⋅ft⁻¹).
"""
@pure pressure(F::FluidState) = F.P

for op ∈ Properties
    @eval @pure $op(F::FluidState) = $op(fluid(F))
end
for op ∈ Intrinsic
    @eval @pure $op(F::FluidState) = $op(temperature(F),fluid(F))
end

@doc """
    viscosity(F::FluidState) = viscosity(temperature(F),fluid(F))

Laminar dynamic vicsosity `μ` at the temperature of `F` (Pa⋅s or slug⋅ft⁻¹⋅s⁻¹).
""" viscosity(F::FluidState)

@doc """
    conductivity(F::FluidState) = conductivity(temperature(F),fluid(F))

Laminar thermal conductivity `k` at the temperature `F` (W⋅m⁻¹⋅K⁻¹ or lb⋅s⁻¹⋅°R⁻¹).
""" conductivity(F::FluidState)

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
    energy(F::FluidState) = energy(temperature(F),fluid(F))

Specific energy `e` at the temperature of `F` (J⋅kg⁻¹ or ft⋅lb⋅slug⁻¹).
""" energy(F::FluidState)

@doc """
    enthalpy(F::FluidState) = enthalpy(temperature(F),fluid(F))

Specific enthalpy `h` at the temperature of `F` (J⋅kg⁻¹ or ft⋅lb⋅slug⁻¹).
""" enthalpy(F::FluidState)

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
@pure density(F::FluidState) = (pressure(F)/temperature(F))/gasconstant(F)

"""
    volume(F::FluidState) = 1/density(F)

Specific volume per mass `v` at a pressure and temperature (m³⋅kg⁻¹ or ft³⋅slug⁻¹).
"""
@pure volume(F::FluidState) = inv(density(F))

"""
    kinematic(F::FluidState) = viscosity(F)/density(F)

Kinematic viscosity ratio `ν` at a pressure and temperature (m²⋅s⁻¹ or ft²⋅s⁻¹).
"""
@pure kinematic(F::FluidState) = viscosity(F)/density(F)

"""
    heatcapacity(F::FluidState) = heatpressure(F)*density(F)

Specific heat per mass at a pressure and temperature (J⋅m⁻³⋅K⁻¹ or lb⋅ft⁻²⋅°R⁻¹).
"""
@pure heatcapacity(F::FluidState) = heatpressure(F)*density(F)

"""
    diffusivity(F::FluidState) = conductivity(F)/heatcapacity(F)

Thermal diffusivity `α` at a pressure and temperature (m²⋅s⁻¹ or ft²⋅s⁻¹).
"""
@pure diffusivity(F::FluidState) = conductivity(F)/heatcapacity(F)

"""
    elasticity(F::FluidState) = heatratio(F)*pressure(F)

Bulk modulus of elasticity `B` at a pressure and temperature (Pa or slug⋅ft⁻¹⋅s⁻²).
"""
@pure elasticity(F::FluidState) = heatratio(F)*pressure(F)

"""
    impedance(F::FluidState) = density(F)*sonicspeed(F)

Specific acoustic resistance at a pressure and temperature (kg⋅m⁻³⋅s⁻¹ or slug⋅ft⁻³⋅s⁻¹).
"""
@pure impedance(F::FluidState) = density(F)*sonicspeed(F)

"""
    intensity(F::FluidState) = pressure(F)^2/impedance(F)

Instantaneous acoustic intensity `I` at a pressure and temperature (W⋅m⁻² or slug⋅s⁻³).
"""
@pure intensity(F::FluidState) = pressure(F)^2/impedance(F)

# global geophysics

"""
    Atmosphere{f,n}

Temperature column of fluid `f` at sea level radius `r` having `n` thermal layers.
"""
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


"""
    Weather{r,g,f,n}

Thermodynamic column state of fluid `f` at sea level radius `r` and gravitational acceleration `g` having `n` thermal `Atmosphere` layers.
Induces derived values `fluid`, `temperature`, `pressure`, `density`, `volume`, `kinematic`, `heatcapacity`, `diffusivity`, `elasticity`, `impedance`, `intensity`, `weight`, `potential`, and inherited values associated with `f::AbstractMole` derivations.
"""
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
(A::Atmosphere)(T,p=atm,r=6.356766e6,g=9.80665) = Weather{r,g}(A,fluid(A)(T,p))
(W::Weather)(hG::Real=0) = W(hG,layer(hG,W))
function (W::Weather)(hG::Real,i)
    h = altgeopotent(hG,W)
    T = _temperature(h,i,W)
    FluidState{fluid(W)}(T,pressure(h,T,i,W))
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

"""
    radius(::Weather)

Sea level radius `r` to planet's center of gravity at `Weather` column location (m or ft).
"""
@pure radius(::Weather{r}=Standard) where r = r

"""
    gravity(::Weather)

Sea level gravitational acceleration `g` at `Weather` column location (m⋅s⁻² or ft⋅s⁻²).
"""
@pure gravity(::Weather{r,g}=Standard) where {r,g} = g

for op ∈ Properties
    @eval @pure $op(W::Weather=Standard) = $op(fluid(W))
end

# hG = geopotential altitude, h = geometric altitude

"""
    altabs(h::Real,W::Weather=Standard) = radius(W)+h

Absolute altitude from planet's center of gravity (m or ft).
"""
@pure @inline altabs(h::Real=0,W::Weather=Standard) = radius(W)+h

"""
    altgeopotent(h::Real,W::Weather=Standard) = h*radius(W)/altabs(h,W)

Geopotential altitude `hG` conversion from geometric altitude (m or ft).
"""
@pure altgeopotent(h::Real,W::Weather=Standard) = (h/altabs(h,W))radius(W)

"""
    altgeometric(hG::Real,W::Weather=Standard) = radius(W)/(radius(W)/hG-1)

Geometric altitude `h` conversion from geopotential altitude (m or ft).
"""
@pure altgeometric(hG::Real,W::Weather=Standard) = (r=radius(W); r/(r/hG-1))

"""
    gravity(h::Real=0,W::Weather=Standard) = gravity(W)*radius(W)^2/altabs(h,W)^2

Gravitational acceleration `g` at altitude `h` of `Weather` column (m⋅s⁻² or ft⋅s⁻²).
"""
@pure gravity(h::Real,W::Weather=Standard) = (gravity(W)*radius(W)^2)/altabs(h,W)^2

"""
    temperature(h::Real=0,::Weather=Standard)

Absolute temperature `T` at geometric altitude `h` of `Weather` location (K or °R).
"""
@pure temperature(h::Real,i,W::Weather=Standard) = _temperature(altgeopotent(h,W),i,W)
@pure function _temperature(hG::Real,i,W::Weather=Standard)
    T0,a0,h0 = W[i]
    return T0+a0*(hG-h0)
end

# k, μ, cᵥ, cₚ, γ, Pr, a, e, h

for op ∈ Intrinsic
    @eval @pure $op(h::Real,i,W::Weather=Standard) = $op(temperature(h,i,W),fluid(W))
end

@doc """
    viscosity(h::Real=0,W::Weather) = viscosity(temperature(h,W),fluid(W))

Laminar dynamic vicsosity `μ` at altitude `h` of `Weather` column (Pa⋅s or slug⋅ft⁻¹⋅s⁻¹).
""" viscosity(h::Real,W::Weather)

@doc """
    conductivity(h::Real=0,W::Weather) = conductivity(temperature(h,W),fluid(W))

Laminar thermal conductivity `k` at altitude `h` of `Weather` (W⋅m⁻¹⋅K⁻¹ or lb⋅s⁻¹⋅°R⁻¹).
""" conductivity(h::Real,W::Weather)

@doc """
    heatvolume(h::Real=0,W::Weather) = heatvolume(temperature(h,W),fluid(W))

Specific heat `cᵥ` at altitude `h` of `Weather` column (J⋅kg⁻¹⋅K⁻¹ or ft⋅lb⋅slug⁻¹⋅°R⁻¹).
""" heatvolume(h::Real,W::Weather)

@doc """
    heatpressure(h::Real=0,W::Weather) = heatpressure(temperature(h,W),fluid(W))

Specific heat `cₚ` at atltitude `h` of `Weather` column (J⋅kg⁻¹⋅K⁻¹ or ft⋅lb⋅slug⁻¹⋅°R⁻¹).
""" heatpressure(h::Real,W::Weather)

@doc """
    heatratio(h::Real=0,W::Weather) = heatratio(temperature(h,W),fluid(W))

Specific heat ratio `γ` at altitude `h` of `Weather` location (dimensionless).
""" heatratio(h::Real,W::Weather)

@doc """
    energy(h::Real=0,W::Weather) = energy(temperature(h,W),fluid(W))

Specific energy `e` at altitude `h` of `Weather` location (J⋅kg⁻¹ or ft⋅lb⋅slug⁻¹).
""" energy(h::Real,W::Weather)

@doc """
    enthalpy(h::Real=0,W::Weather) = enthalpy(temperature(h,W),fluid(W))

Specific enthalpy `h` at altitude of `Weather` location (J⋅kg⁻¹ or ft⋅lb⋅slug⁻¹).
""" enthalpy(h::Real,W::Weather)

@doc """
    freedom(h::Real=0,W::Weather) = freedom(temperature(h,W),fluid(W))

Degrees of freedom `f` at altitude `h` of `Weather location (dimensionless).
""" freedom(h::Real,W::Weather)

@doc """
    prandtl(h::Real=0,W::Weather) = prandtl(temperature(h,W),fluid(W))

Prandtl number at altitude `h` of `Weather` location (dimensionless).
""" prandtl(h::Real,W::Weather)

@doc """
    sonicspeed(h::Real=0,W::Weather) = sonicspeed(temperature(h,W),fluid(W))

Speed of sound wave disturbance at altitude `h` of `Weather` location (m⋅s⁻¹ or ft⋅s⁻¹).
""" sonicspeed(h::Real,W::Weather)

"""
    pressure(h::Real=0,W::Weather=Standard) = pressure(W(h))

Absolute force per unit area `P`  at altitude `h` of `Weather` column (Pa or slug⋅ft⁻¹⋅s⁻²).
"""
@pure function pressure(h::Real,i,W::Weather=Standard)
    hG = altgeopotent(h,W)
    pressure(h,_temperature(hG,i,W),i,W)
end
@pure function pressure(hG::Real,T,i,W::Weather=Standard)
    g,R = gravity(W),gasconstant(W)
    T0,a,h0,p = W[i]
    p*(if iszero(a)
        exp((-g/R)*(hG-h0)/T)
    else
        (T/T0)^((-g/R)/a)
    end)
end

"""
    density(h::Real=0,W::Weather=Standard) = density(W(h))

Inertial mass per volume `ρ` at altitude `h` of `Weather` location (kg⋅m⁻³ or slugs⋅ft⁻³).
"""
@pure function density(h::Real,i,W::Weather=Standard)
    hG = altgeopotent(h,W)
    density(h,_temperature(hG,i,W),i,W)
end
@pure function density(hG::Real,T,i,W::Weather=Standard)
    g,R = gravity(W),gasconstant(W)
    T0,a,h0,_,ρ = W[i]
    ρ*(if iszero(a)
        exp((-g/R)*(hG-h0)/T)
    else
        (T/T0)^((-g/R)/a-1)
    end)
end

"""
    kinematic(h::Real=0,W::Weather=Standard) = kinematic(W(h))

Kinematic viscosity ratio `ν` at altitude `h` of `Weather` location (m²⋅s⁻¹ or ft²⋅s⁻¹).
"""
@pure function kinematic(h::Real,i,W::Weather=Standard)
    hG = altgeopotent(h,W)
    kinematic(h,_temperature(hG,i,W),i,W)
end
@pure kinematic(hG::Real,T,i,W::Weather=Standard) = viscosity(T,fluid(W))/density(hG,T,i,W)

"""
    heatcapacity(h::Real=0,W::Weather=Standard) = heatcapacity(W(h))

Specific heat per mass at altitude `h` of `Weather` location (J⋅m⁻³⋅K⁻¹ or lb⋅ft⁻²⋅°R⁻¹).
"""
@pure function heatcapacity(h::Real,i,W::Weather=Standard)
    hG = altgeopotent(h,W)
    T = _temperature(hG,i,W)
    heatpressure(T,fluid(W))*density(hG,T,i,W)
end

"""
    diffusivity(h::Real=0,W::Weather=Standard) = diffusivity(W(h))

Thermal diffusivity `α` at altitude `h` of `Weather` location (m²⋅s⁻¹ or ft²⋅s⁻¹).
"""
@pure function diffusivity(h::Real,i,W::Weather=Standard)
    F = fluid(W)
    hG = altgeopotent(h,W)
    T = _temperature(h,i,W)
    conductivity(T,F)/heatpressure(T,F)/density(hG,T,i,W)
end

"""
    elasticity(h::Real=0,W::Weather=Standard) = elasticity(W(h))

Bulk modulus of elasticity `B` at altitude `h` of `Weather` location (Pa or slug⋅ft⁻¹⋅s⁻²).
"""
@pure function elasticity(h,i,W::Weather=Standard)
    hG = altgeopotent(h,W)
    T = _temperature(hG,i,W)
    heatratio(T,fluid(W))*pressure(hG,T,i,W)
end

"""
    impedance(h::Real=0,W::Weather=Standard) = impedance(W(h))

Specific acoustic resistance at altitude `h` of `Weather` (kg⋅m⁻³⋅s⁻¹ or slug⋅ft⁻³⋅s⁻¹).
"""
@pure function impedance(h::Real,i,W::Weather=Standard)
    hG = altgeopotent(h,W)
    T = _temperature(hG,i,W)
    density(hG,T,i,W)*sonicspeed(T,fluid(W))
end

"""
    intensity(h::Real=0,W::Weather=Standard) = intensity(W(h))

Instantaneous intensity `I` at altitude `h` of `Weather` at location (W⋅m⁻² or slug⋅s⁻³).
"""
@pure function intensity(h::Real,i,W::Weather=Standard)
    hG = altgeopotent(h,W)
    T = _temperature(hG,i,W)
    g,R = gravity(W),gasconstant(W)
    T0,a,h0,p,ρ = W[i]
    (p^2/ρ)*(if iszero(a)
        exp((-g/R)*(hG-h0)/T)
    else
        t,gRa = T/T0,(-g/R)/a
        t^2gRa/t^(gRa-1)
    end)/sonicspeed(T,fluid(W))
end

# Grashof number

@pure function grashof(h::Real,i,W::Weather=Standard)
    hG = altgeopotent(h,W)
    T = _temperature(hG,i,W)
    gravity(h,W)*(temperature(W)-T)*(h^3)/(T*kinematic(hG,T,i,W)^2)
end

"""
    weight(h::Real=0,W::Weather=Standard) = density(h,W)*gravity(h,W)

Specific weight at altitude `h` of `Weather` location (kg⋅m⁻²⋅s⁻² or slugs⋅ft⁻²⋅s⁻²).
"""
@pure weight(h::Real,i,W::Weather=Standard) = density(h,i,W)*gravity(h,W)

"""
    volume(h::Real=0,W::Weather=Standard) = volume(W(h))

Specific volume per mass `v` at altitude `h` of `Weather` location (m³⋅kg⁻¹ or ft³⋅slug⁻¹).
"""
@pure volume(h::Real,i,W::Weather=Standard) = inv(density(h,i,W))

"""
    potential(h::Real=0,W::Weather=Standard) = gravity(h,W)*h

Specifc gravitational potential energy `g` at altitude `h` of `Weather` (m⋅s⁻² or ft⋅s⁻²).
"""
@pure potential(h::Real,i,W::Weather=Standard) = gravity(h,W)*h

# common interface

for op ∈ (:temperature,:pressure,:density,:weight,:volume,:potential,:impedance,:grashof,:diffusivity,:intensity,:heatcapacity,:kinematic,:elasticity,Intrinsic...)
    opratio = Symbol(op,:ratio)
    @eval begin
        export $op
        @pure $op(W::Weather=Standard) = $op(0,W)
        @pure $op(h::Real,W::Weather=Standard) = $op(h,layer(h,W),W)
    end
    op ∉ (:heatcapacity,:grashof,:potential) && @eval begin
        export $opratio
        @pure $opratio(h::Real,W::Weather=Standard) = $opratio(h,layer(h,W),W)
        @pure $opratio(h::Real,i,W::Weather=Standard) = $op(h,i,W)/$op(W)
    end
end

# data

include("planets.jl")

end # module
