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

using UnitSystems
export Metric, English
import UnitSystems: Constants, molarmass, units
for op ∈ (Constants...,UnitSystems.Physics...,UnitSystems.Convert...)
    @eval import UnitSystems.$op
end

const Properties = (:molecularmass,:gasconstant,UnitSystems.Constants...)
const Intrinsic = (:viscosity,:thermalconductivity,:heatvolume,:heatpressure,:heatratio,:prandtl,:sonicspeed,:freedom,:specificenergy,:specificenthalpy)

include("chemistry.jl")

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

# global geophysics

"""
    Atmosphere{n,U}

Temperature column of planet in unit system `U` having `n` thermal layers.
"""
struct Atmosphere{n,U}
    a::Values{n,Float64} # lapse rate
    h::Values{n,Float64} # altitude
    m::Values{n,Float64} # molar rate
    Atmosphere{U}(a::Values{n},h::Values{n},m::Values{n}) where {n,U} = new{n,U}(a,h,m)
end

@pure units(::Atmosphere{n,U}) where {n,U} = U
Atmosphere(a::Values{n},h::Values{n}) where n = Atmosphere{Metric}(a,h)
Atmosphere{U}(a::Values{n},h::Values{n}) where {n,U} = Atmosphere{U}(a,h,zeros(Values{n}))
(U::UnitSystem)(A::Atmosphere{n,S}) where {n,S} = Atmosphere{U}(lapserate.(A.a,Ref(U),Ref(S)),length.(A.h,Ref(U),Ref(S)))

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
Induces derived values `fluid`, `temperature`, `pressure`, `density`, `specificvolume`, `kinematic`, `heatcapacity`, `thermaldiffusivity`, `elasticity`, `specificimpedance`, `intensity`, `specificweight`, `geopotential`, and inherited values associated with `f::AbstractMole` derivations.
"""
struct Weather{r,g,f,n,U} # radius, acceleration
    A::Atmosphere{n,U} # altitude lapse rate
    T::Values{n,Float64} # temperature
    p::Values{n,Float64} # pressure
    ρ::Values{n,Float64} # density
    Weather{r,g,f}(A::Atmosphere{n,U},T::Values{n},p::Values{n},ρ::Values{n}) where {r,g,f,n,U} = new{r,g,f,n,U}(A,T,p,ρ)
end

function Weather{r,g}(A::Atmosphere{n,U},F::FluidState{f}) where {r,g,f,n,U}
    T = zeros(Variables{n,Float64})
    p = zeros(Variables{n,Float64})
    ρ = zeros(Variables{n,Float64})
    T0,p0 = temperature(F,U),pressure(F,U)
    R,γ = gasconstant(F,U),heatratio(F,U)
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
    Weather{r,g,f}(A,Values(T),Values(p),Values(ρ))
end

(U::UnitSystem)(W::Weather{r,g,f,n,S}) where {r,g,f,n,S} = Weather{r,g,f}(U(W.A),temperature.(W.T,Ref(U),Ref(S)),pressure.(W.p,Ref(U),Ref(S)),density.(W.ρ,Ref(U),Ref(S)))
(A::Atmosphere)(F::FluidState=Air(288.16,atm,US(A)),r=6.356766e6,g=g₀) = Weather{r,g}(A,F)
(A::Atmosphere)(T,p=atm,r=6.356766e6,g=g₀) = Weather{r,g}(A,Air(T,p,units(A)))
(W::Weather)(hG::Real=0) = W(hG,layer(hG,W))
function (W::Weather)(hG::Real,i)
    h = altgeopotent(hG,W)
    T = _temperature(h,i,W)
    FluidState{fluid(W),units(W)}(T,pressure(h,T,i,W))
end

function display(W::Weather)
    println(typeof(W))
    println(" a = ",W.A.a)
    println(" h = ",W.A.h)
    println(" T = ",W.T)
    println(" P = ",W.p)
    println(" ρ = ",W.ρ)
end

@pure @inline function Base.getindex(W::Weather,i::Int,U::UnitSystem=units(W))
    T,a,h,p,ρ = W.T[i],W.A.a[i],W.A.h[i],W.p[i],W.ρ[i]; S = units(W)
    temperature(T,U,S),lapserate(a,U,S),length(h,U,S),pressure(p,U,S),density(ρ,U,S)
end
@pure @inline Base.getindex(W::Weather,::Val{i}) where i = getindex(W,i)
@pure @inline lapserate(h::Real,W::Weather=Standard) = W.A.a[layer(h,W)]
@pure @inline layer(h::Real,W::Weather=Standard) = h≤W.A.h[1] ? 1 : (i=findfirst(x->x≥h,W.A.h); isnothing(i) ? length(W.A.h) : i-1)
@pure @inline layer(h::Real,W::Weather,U::US) = layer(length(h,units(W),U),W)

@pure units(::Weather{r,g,f,n,U}) where {r,g,f,n,U} = U
@pure fluid(::Weather{r,g,f}=Standard) where {r,g,f} = f

"""
    radius(::Weather)

Sea level radius `r` to planet's center of gravity at `Weather` column location (m or ft).
"""
@pure radius(W::Weather{r}=Standard,U::US=units(W)) where r = length(r,U,units(W))

"""
    gravity(::Weather)

Sea level gravitational acceleration `g` at `Weather` column location (m⋅s⁻² or ft⋅s⁻²).
"""
@pure gravity(W::Weather{r,g}=Standard,U::US=units(W)) where {r,g} = acceleration(g,U,units(W))

@pure molecularmass(W::Weather=Standard,U::US=units(W)) = molecularmass(fluid(W),U)
@pure gasconstant(W::Weather=Standard,U::US=units(W)) = gasconstant(fluid(W),U)

# hG = geopotential altitude, h = geometric altitude

"""
    altabs(h::Real,W::Weather=Standard) = radius(W)+h

Absolute altitude from planet's center of gravity (m or ft).
"""
@pure @inline altabs(h::Real=0,W::Weather=Standard,U::US=units(W)) = radius(W,U)+h
@pure @inline altabs(h::Real,W::Weather,U::US,S::US) = altabs(length(h,U,S),W,U)

"""
    altgeopotent(h::Real,W::Weather=Standard) = h*radius(W)/altabs(h,W)

Geopotential altitude `hG` conversion from geometric altitude (m or ft).
"""
@pure altgeopotent(h::Real,W::Weather=Standard,U::US=US(W)) = (h/altabs(h,W,U))radius(W,U)
@pure altgeopotent(h::Real,W::Weather,U::US,S::US) = altgeopotent(length(h,U,S),W,U)

"""
    altgeometric(hG::Real,W::Weather=Standard) = radius(W)/(radius(W)/hG-1)

Geometric altitude `h` conversion from geopotential altitude (m or ft).
"""
@pure altgeometric(hG::Real,W::Weather=Standard,U::US=US(W)) = (r=radius(W,U); r/(r/hG-1))
@pure altgeometric(hG::Real,W::Weather,U::US,S::US) = altgeometric(length(hG,U,S),W,U)

"""
    gravity(h::Real=0,W::Weather=Standard) = gravity(W)*radius(W)^2/altabs(h,W)^2

Gravitational acceleration `g` at altitude `h` of `Weather` column (m⋅s⁻² or ft⋅s⁻²).
"""
@pure gravity(h::Real,W::Weather=Standard,U::US=US(W)) = (gravity(W,U)*radius(W,U)^2)/altabs(h,W,U)^2
@pure gravity(h::Real,W::Weather,U::US,S::US) = gravity(length(h,U,S),W,U)

"""
    temperature(h::Real=0,::Weather=Standard)

Absolute temperature `T` at geometric altitude `h` of `Weather` location (K or °R).
"""
@pure temperature(h::Real,i,W::Weather=Standard,U::US=US(W)) = _temperature(altgeopotent(h,W,U),i,W,U)
@pure function _temperature(hG::Real,i,W::Weather=Standard,U::US=units(W))
    T0,a0,h0 = W[i,U]
    return T0+a0*(hG-h0)
end

# k, μ, cᵥ, cₚ, γ, Pr, a, e, h

for op ∈ Intrinsic
    @eval @pure $op(h::Real,i,W::Weather=Standard,U::US=US(W)) = $op(temperature(h,i,W,U),fluid(W),U)
    @eval @pure $op(h::Real,i,W::Weather,U::US,S::US) = $op(length(h,U,S),i,W,U)
end

@doc """
    viscosity(h::Real=0,W::Weather) = viscosity(temperature(h,W),fluid(W))

Laminar dynamic vicsosity `μ` at altitude `h` of `Weather` column (Pa⋅s or slug⋅ft⁻¹⋅s⁻¹).
""" viscosity(h::Real,W::Weather)

@doc """
    thermalconductivity(h::Real=0,W::Weather) = thermalconductivity(temperature(h,W),fluid(W))

Laminar thermal conductivity `k` at altitude `h` of `Weather` (W⋅m⁻¹⋅K⁻¹ or lb⋅s⁻¹⋅°R⁻¹).
""" thermalconductivity(h::Real,W::Weather)

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
    specificenergy(h::Real=0,W::Weather) = specificenergy(temperature(h,W),fluid(W))

Specific energy `e` at altitude `h` of `Weather` location (J⋅kg⁻¹ or ft⋅lb⋅slug⁻¹).
""" specificenergy(h::Real,W::Weather)

@doc """
    specificenthalpy(h::Real=0,W::Weather) = specificenthalpy(temperature(h,W),fluid(W))

Specific enthalpy `h` at altitude of `Weather` location (J⋅kg⁻¹ or ft⋅lb⋅slug⁻¹).
""" specificenthalpy(h::Real,W::Weather)

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
@pure function pressure(h::Real,i,W::Weather=Standard,U::US=units(W))
    hG = altgeopotent(h,W,U)
    pressure(h,_temperature(hG,i,W,U),i,W,U)
end
@pure function pressure(hG::Real,T,i,W::Weather=Standard,U::US=units(W))
    g,R = gravity(W,U),gasconstant(W,U)
    T0,a,h0,p = W[i,U]
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
@pure function density(h::Real,i,W::Weather=Standard,U::US=units(W))
    hG = altgeopotent(h,W,U)
    density(h,_temperature(hG,i,W,U),i,W,U)
end
@pure density(hG::Real,T,i,W::Weather,U::US,S::US) = density(length(h,U,S),temperature(T,U,S),i,W,U)
@pure function density(hG::Real,T,i,W::Weather=Standard,U::US=units(W))
    g,R = gravity(W,U),gasconstant(W,U)
    T0,a,h0,_,ρ = W[i,U]
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
@pure function kinematic(h::Real,i,W::Weather=Standard,U::US=units(W))
    hG = altgeopotent(h,W,U)
    kinematic(h,_temperature(hG,i,W,U),i,W,U)
end
@pure kinematic(hG::Real,T,i,W::Weather,U::US,S::US) = kinematic(length(h,U,S),temperature(T,U,S),i,W,U)
@pure kinematic(hG::Real,T,i,W::Weather=Standard,U::US=units(W)) = viscosity(T,fluid(W),U)/density(hG,T,i,W,U)

"""
    heatcapacity(h::Real=0,W::Weather=Standard) = heatcapacity(W(h))

Specific heat per mass at altitude `h` of `Weather` location (J⋅m⁻³⋅K⁻¹ or lb⋅ft⁻²⋅°R⁻¹).
"""
@pure function heatcapacity(h::Real,i,W::Weather=Standard,U::US=units(W))
    hG = altgeopotent(h,W,U)
    T = _temperature(hG,i,W,U)
    heatpressure(T,fluid(W),U)*density(hG,T,i,W,U)
end

"""
    thermaldiffusivity(h::Real=0,W::Weather=Standard) = thermaldiffusivity(W(h))

Thermal diffusivity `α` at altitude `h` of `Weather` location (m²⋅s⁻¹ or ft²⋅s⁻¹).
"""
@pure function thermaldiffusivity(h::Real,i,W::Weather=Standard,U::US=units(W))
    F,hG = fluid(W),altgeopotent(h,W,U)
    T = _temperature(hG,i,W,U)
    thermalconductivity(T,F,U)/heatpressure(T,F,U)/density(hG,T,i,W,U)
end

"""
    elasticity(h::Real=0,W::Weather=Standard) = elasticity(W(h))

Bulk modulus of elasticity `B` at altitude `h` of `Weather` location (Pa or slug⋅ft⁻¹⋅s⁻²).
"""
@pure function elasticity(h,i,W::Weather=Standard,U::US=units(W))
    hG = altgeopotent(h,W,U)
    T = _temperature(hG,i,W,U)
    heatratio(T,fluid(W),U)*pressure(hG,T,i,W,U)
end

"""
    specificimpedance(h::Real=0,W::Weather=Standard) = impedance(W(h))

Specific acoustic resistance at altitude `h` of `Weather` (kg⋅m⁻³⋅s⁻¹ or slug⋅ft⁻³⋅s⁻¹).
"""
@pure function specificimpedance(h::Real,i,W::Weather=Standard,U::US=units(W))
    hG = altgeopotent(h,W,U)
    T = _temperature(hG,i,W,U)
    density(hG,T,i,W,U)*sonicspeed(T,fluid(W),U)
end

"""
    intensity(h::Real=0,W::Weather=Standard) = intensity(W(h))

Instantaneous intensity `I` at altitude `h` of `Weather` at location (W⋅m⁻² or slug⋅s⁻³).
"""
@pure function intensity(h::Real,i,W::Weather=Standard,U::US=units(W))
    hG = altgeopotent(h,W,U)
    T = _temperature(hG,i,W,U)
    g,R = gravity(W,U),gasconstant(W,U)
    T0,a,h0,p,ρ = W[i,U]
    (p^2/ρ)*(if iszero(a)
        exp((-g/R)*(hG-h0)/T)
    else
        t,gRa = T/T0,(-g/R)/a
        t^2gRa/t^(gRa-1)
    end)/sonicspeed(T,fluid(W),U)
end

# Grashof number

@pure function grashof(h::Real,i,W::Weather=Standard,U::US=units(W))
    hG = altgeopotent(h,W,U)
    T = _temperature(hG,i,W,U)
    gravity(h,W,U)*(temperature(W,U)-T)*(h^3)/(T*kinematic(hG,T,i,W,U)^2)
end

"""
    specificweight(h::Real=0,W::Weather=Standard) = density(h,W)*gravity(h,W)

Specific weight at altitude `h` of `Weather` location (kg⋅m⁻²⋅s⁻² or slugs⋅ft⁻²⋅s⁻²).
"""
@pure specificweight(h::Real,i,W::Weather=Standard,U::US=US(W)) = density(h,i,W,U)*gravity(h,W,U)

"""
    specificvolume(h::Real=0,W::Weather=Standard) = specificvolume(W(h))

Specific volume per mass `v` at altitude `h` of `Weather` location (m³⋅kg⁻¹, ft³⋅slug⁻¹).
"""
@pure specificvolume(h::Real,i,W::Weather=Standard,U::US=US(W)) = inv(density(h,i,W,U))

"""
    geopotential(h::Real=0,W::Weather=Standard) = gravity(h,W)*h

Specifc gravitational potential energy `g` at altitude `h` of `Weather` (m²⋅s⁻², ft²⋅s⁻²).
"""
@pure geopotential(h::Real,i,W::Weather=Standard,U::US=US(W)) = gravity(h,W,U)*h

# common interface

for op ∈ (:temperature,:pressure,:density,:specificweight,:specificvolume,:geopotential,:specificimpedance,:grashof,:thermaldiffusivity,:intensity,:heatcapacity,:kinematic,:elasticity,Intrinsic...)
    opratio = Symbol(op,:ratio)
    @eval begin
        export $op
        @pure $op(W::Weather=Standard,U::US=Metric) = $op(0,W,U)
        @pure $op(h::Real,W::Weather=Standard,U::US=US(W)) = $op(h,layer(h,W,U),W,U)
        @pure $op(h::Real,W::Weather,U::US,S::US) = $op(length(h,U,S),W,U)
        @pure $op(h,i,W::Weather,U::US,S::US) = $op(length(h,U,S),i,W,U)
    end
    op ∉ (:heatcapacity,:grashof,:geopotential) && @eval begin
        export $opratio
        @pure $opratio(h::Real,W::Weather=Standard,U::US=US(W)) = $opratio(h,layer(h,W,U),W,U)
        @pure $opratio(h::Real,i,W::Weather=Standard,U::US=US(W)) = $op(h,i,W,U)/$op(W,U)
        @pure $opratio(h::Real,W::Weather,U::US,S::US) = $opratio(length(h,U,S),W,U)
        @pure $opratio(h::Real,i,W::Weather,U::US,S::US) = $opratio(length(h,U,S),i,W,U)
    end
end

# data

include("planets.jl")

end # module
