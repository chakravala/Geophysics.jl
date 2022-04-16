module Geophysics

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
# (                            )
# )\ )      (               ( /(  (        (
#(()/(     ))\  (    `  )   )\()) )\ )  (  )\   (   (
# /(_))_  /((_) )\   /(/(  ((_)\ (()/(  )\((_)  )\  )\
#(_)) __|(_))  ((_) ((_)_\ | |(_) )(_))((_)(_) ((_)((_)
#  | (_ |/ -_)/ _ \ | '_ \)| ' \ | || |(_-<| |/ _| (_-<
#   \___|\___|\___/ | .__/ |_||_| \_, |/__/|_|\__| /__/
#                   |_|             |_|

using AbstractTensors, LinearAlgebra
import Base: @pure, show, display

export FluidState, Atmosphere, Weather, Planet

export fluid, radius, gravity, layer

#decibel(I₁,I₂=intensity()) = 10log10(I₁/I₂)
gage(P::Real,P0::Real=pressure()) = P-P0

export Metric, English, British
import UnitSystems
import UnitSystems: units
const usingSimilitude = UnitSystems.similitude()

if usingSimilitude
import UnitSystems
using Similitude
import Similitude: Constants, molarmass
for op ∈ (Constants...,UnitSystems.Physics...,UnitSystems.Convert...)
    @eval import Similitude.$op
end
else
using UnitSystems
import UnitSystems: Quantity, Constants, molarmass, normal
for op ∈ (Constants...,UnitSystems.Physics...,UnitSystems.Convert...)
    @eval import UnitSystems.$op
end
end

const Properties = (:molecularmass,:gasconstant,UnitSystems.Constants...)
const Intrinsic = (:viscosity,:thermalconductivity,:heatvolume,:heatpressure,:heatratio,:prandtl,:sonicspeed,:freedom,:specificenergy,:specificenthalpy)

export flattening, semimajor, semiminor, period, gravitation, frequency, angularfrequency
export eccentricity, eccentricity2, lineareccentricity, aspectratio, radius, speed, mass
export latitudegeodetic, latitudegeocentric, latitudeparametric
export deflectiongeodetic, deflectiongeocentric, deflection
export centripetal, oblateness, q0, q01, q, q1, dynamicformfactor, secondzonalharmonic
export gravitycomponents, gravitygeodetic, radiusgeodetic

include("chemistry.jl")

# global geophysics

"""
    Planet{f,a,t,Gm}

Celestial `Planet` ellipsoidal reference body with `flattening` of `f`, `semimajor` radius of `a`, sidereal `period` of `t`, and standard `gravitation` parameter `Gm`.
Further derived parameters include `semiminor`, `frequency`, `angularfrequency`, `eccentricity`, `eccentricity2`, `lineareccentricity`, `aspectratio`, `oblateness`, `dynamicformfactor`, `secondzonalharmonic`, `latitudegeodetic`, `latitudegeocentric`, `latitudeparametric`, `radius`, `radiusgeodetic`, `gravity`, `gravitygeodetic`, and `gravitycomponents`.
"""
struct Planet{f,a,t,Gm} end # WGS 84 spheroid
Planet(f,a,t,Gm,U=Metric) = Planet{f,Quantity(L,U,a),Quantity(T,U,t),Quantity(F/M*L^2,U,Gm)}()
const Earth = Planet(1/298.257223563,6378137.0,86164.098903691,3.986004418e14)

"""
    flattening(P::Planet) = 1-semiminor(P)/semimajor(P)

Oblate `flattening` parameter of a `Planet` reference ellipsoid (dimensionless).

```Julia
julia> 1/flattening(Earth)
$(1/flattening(Earth))
```
"""
@pure flattening(::Planet{f}=Earth) where {f} = f

"""
    semimajor(P::Planet) = semiminor(P)/(1-flattening(P))

Equatorial radius of oblate `Planet` reference ellipsoid (m).

```Julia
julia> semimajor(Earth) # m
$(semimajor(Earth))
```
"""
@pure semimajor(P::Planet{f,a}=Earth,U::US=Metric) where {f,a} = a*length(Metric,U)

"""
    period(P::Planet) = 2π/angularfrequency(P)

Sidereal rotational `period` of `Planet` (s).

```Julia
julia> period(Earth) # s
$(period(Earth))
```
"""
@pure period(::Planet{f,a,t}=Earth,U::US=Metric) where {f,a,t} = t*time(Metric,U)

"""
    gravitation(P::Planet,U::UnitSystem) = mass(P)*gravitation(U)

Standard `gravitation` parameter for a celestial `Planet` body (m³⋅s⁻²).

```Julia
julia> gravitation(Earth) # m³⋅s⁻²
$(gravitation(Earth))
```
"""
@pure gravitation(::Planet{f,a,t,Gm}=Earth) where {f,a,t,Gm} = Gm
@pure gravitation(P::Planet,U::US) = gravitation(P)/(length(U,Metric)*specificenergy(U,Metric))

"""
    mass(P::Planet,U::UnitSystem) = gravitation(P,U)/newton(U)

Inertial `mass` of a celestial `Planet` body (kg).

```Julia
julia> mass(Earth) # kg
$(mass(Earth))
```
"""
@pure mass(P::Planet,U::US=Metric) = gravitation(P,U)/gravitation(U)

"""
    frequency(P::Planet) = 1/period(P)

Cyclic `frequency` of sidereal `Planet` rotation (Hz).

```Julia
julia> frequency(Earth) # Hz
$(frequency(Earth))
```
"""
@pure frequency(P::Planet,U::US=Metric) = 1/period(P,U)

"""
    angularfrequency(P::Planet) = 2π/period(P)

Rate of `angularfrequency` of sidereal `Planet` rotation (rad⋅s⁻¹).

```Julia
julia> angularfrequency(Earth) # rad⋅s⁻¹
$(angularfrequency(Earth))
```
"""
@pure angularfrequency(P::Planet,U::US=Metric) = 2π/period(P,U)

@pure meanradius(P::Planet,U::US=Metric) = radius_fast(asin(sqrt(1/3)),P,U)
@pure radius_fast(θ,P::Planet,U::US=Metric) = semimajor(P,U)*(1-flattening(P)*sin(θ)^2)

"""
    semiminor(P::Planet) = semimajor(P)*(1-flattening(P))

Polar radius of oblate `Planet` reference ellipsoid (m).

```Julia
julia> semiminor(Earth) # m
$(semiminor(Earth))
```
"""
@pure semiminor(P::Planet,U::US=Metric) = radius_fast(π/2,P,U)

"""
    eccentricity(P::Planet) = sqrt(flattening(P)*(2-flattening(P)))

Elliptic `eccentricity` of oblate `Planet` reference ellipsoid (dimensionless).

```Julia
julia> eccentricity(Earth)
$(eccentricity(Earth))
```
"""
@pure eccentricity(::Planet{f}) where f = sqrt(f*(2-f))

"""
    eccentricity2(P::Planet) = eccentricity(P)/sqrt(1-eccentricity(P))

Second `eccentricity2` of oblate `Planet reference ellipsoid (dimensionless).

```Julia
julia> eccentricity2(Earth)
$(eccentricity2(Earth))
```
"""
@pure eccentricity2(P::Planet{f}) where f = eccentricity(P)/(1-f)

"""
    lineareccentricity(P::Planet) = semimajor(P)*eccentricity(P)

Distance from center to foci of oblate `Planet` reference ellipsoid (m).
"""
@pure lineareccentricity(P::Planet,U::US=Metric) = semimajor(P,U)*eccentricity(P)

"""
    aspectratio(P::Planet) = semiminor(P)/semimajor(P)

Value of `aspectratio` or `semiminor` per `semimajor` of `Planet` (dimensionless).

```Julia
julia> aspectratio(Earth)
$(aspectratio(Earth))
```
"""
@pure aspectratio(P::Planet) = semiminor(P)/semimajor(P)

@pure authalicradius(P::Planet,U::US=Metric) = sqrt((semimajor(P,U)^2+semiminor(P,U)^2*atanh(eccentricity(P))/eccentricity(P))/2)

"""
    latitudegeodetic(θ,P::Planet) = atan(tan(θ)/(1-flattening(P))^2)

Convert geocentric latitude `θ` to geodetic latitude on `Planet` surface (rad).
"""
@pure latitudegeodetic(θ,P::Planet=Earth) = atan(tan(θ)/(1-flattening(P))^2)
@pure deflectiongeodetic(θ,P=Earth) = latitudegeodetic(θ,P)-θ

"""
    latitudegeocentric(ϕ,P::Planet) = atan(tan(ϕ)*(1-flattening(P))^2)

Convert geodetic latitude `ϕ` to geocentric latitude on `Planet` surface (rad).
"""
@pure latitudegeocentric(ϕ,P::Planet=Earth) = atan(tan(ϕ)*(1-flattening(P))^2)
@pure deflectiongeocentric(ϕ,P=Earth) = latitudegeocentric(ϕ,P)-ϕ

"""
    latitudeparametric(ϕ,P::Planet) = atan(tan(ϕ)*(1-flattening(P)))

Convert geodetic latitude `ϕ` to parametric latitude on surrounding sphere (rad).
"""
@pure latitudeparametric(ϕ,P::Planet=Earth) = atan(tan(ϕ)*(1-flattening(P)))

"""
    deflection(h,ϕ,P::Planet)

Approximation for `deflection` angle between geodetic and geocentric latitudes (rad).
"""
@pure function deflection(h,ϕ,P::Planet=Earth,U::US=Metric)
    f = flattening(P)
    f*sin(2ϕ)*(1-f/2-h/radiusgeodetic(ϕ,P,U))
end # geodetic: h, ϕ

"""
    latitudegeocentric(h,ϕ,P::Planet) = ϕ-deflection(h,ϕ,P)

Convert geodetic latitude `ϕ` to geocentric latitude at altitude `h` (rad).
"""
@pure latitudegeocentric(h,ϕ,P::Planet=Earth,U::US=Metric) = ϕ-deflection(h,ϕ,P,U)

"""
    radius(θ,P::Planet) = 1/sqrt((cos(θ)/semimajor(P,U))^2+(sin(θ)/semiminor(P,U))^2)

Parametrized `radius` of `Planet` in terms of geocentric latitude `θ` coordinate (m).
"""
@pure radius(θ,P::Planet,U::US=Metric) = 1/sqrt((cos(θ)/semimajor(P,U))^2+(sin(θ)/semiminor(P,U))^2)

"""
    radiusgeodetic(ϕ,P::Planet)

Approximated `radius` of `Planet` in terms of geodetic latitude `ϕ` coordinate (m).
"""
@pure function radiusgeodetic(ϕ,P::Planet=Earth,U::US=Metric)
    f,a = flattening(P),semimajor(P,U)
    a*(1-f/2*(1-cos(2ϕ))+5f^2/16*(1-cos(4ϕ)))
end # approximation, geodetic: λd

@pure _speed(θ,P::Planet,U::US=Metric) = radius(θ,P,U)*angularfrequency(P,U)
"""
    speed(θ,P::Planet) = cos(θ)*radius(θ,P)*angularfrequency(P)

Velocity component due to sidereal rotation of `Planet` at latitude `θ` (m⋅s⁻¹).

```Julia
julia> speed(0,Earth)
$(speed(0,Earth))
```
"""
@pure speed(θ,P::Planet,U::US=Metric) = _speed(θ,P,U)*cos(θ)

@pure _centripetal(θ,P::Planet=Earth,U::US=Metric) = _speed(θ,P,U)*angularfrequency(P,U)
"""
    centripetal(θ,P::Planet) = speed(θ,P)*angularfrequency(P)

Acceleration due to `centripetal` force on `Planet` at geocentric latitude `θ` (m⋅s⁻²).
"""
@pure centripetal(θ,P::Planet=Earth,U::US=Metric) = _centripetal(θ,P,U)*cos(θ)

"""
    gravity(P::Planet) = gravitation(P,U)/semimajor(P,U)^2

Spherical estimate of gravitational acceleration at equator of `Planet` (m⋅s⁻²)
"""
@pure gravity(P::Planet,U::US=Metric) = gravitation(P,U)/semimajor(P,U)^2

"""
    oblateness(θ,P::Planet) = radius(θ,P)*angularfrequency(P)^2/gravity(P)

Constant of `oblateness` at latitude `θ` for gravitating `Planet` ellipse.
"""
@pure oblateness(θ,P::Planet=Earth,U::US=Metric) = _centripetal(θ,P,U)/gravity(P,U)/gravity(U)

"""
    oblateness(P::Planet) = radius(π/2,P)*angularfrequency(P)^2/gravity(P)

Constant of `oblateness` at latitude `π/2` for gravitating `Planet` ellipse.

```Julia
julia> oblateness(Earth) # θ ≡ π/2
$(oblateness(Earth))
```
""" # normal gravity constant
@pure oblateness(P::Planet=Earth,U::US=Metric) = oblateness(π/2,P,U)

# Hirvonen zonal notation
@pure q0(P::Planet) = (e2=eccentricity2(P); ((1+3/e2^2)*atan(e2)-3/e2)/2)
@pure q01(P::Planet) = (e2=eccentricity2(P); 3((1+1/e2^2)*(1-atan(e2)/e2))-1)
@pure q(u,P,U) = (E=lineareccentricity(P,U)/u; ((1+3/E^2)*atan(E)-3/E)/2)
@pure q1(u,P,U) = (E=lineareccentricity(P,U)/u; 3((1+E^-2)*(1-atan(E)/E))-1)

# should be zero maybe, but due to 0/0 are set to 1
@pure q0(P::Planet{0}) = 1 # ?
@pure q01(P::Planet{0}) = 1 # ?
@pure q(u,P::Planet{0},U) = 1 # ?
@pure q1(u,P::Planet{0},U) = 1 # ?

"""
    dynamicormfactor(::Planet)

Celestial body's second `dynamicformfactor` in nodal precession (dimensionless).

```Julia
julia> dynamicformfactor(Earth)
$(dynamicformfactor(Earth))
```
"""
@pure dynamicformfactor(P::Planet{f}=Earth) where f = (1-2oblateness(P)*eccentricity2(P)/15q0(P))*f*(2-f)/3
@pure dynamicformfactor(P::Planet{0}) = 0

"""
    secondzonhalharmonic(P::Planet) = -dynamicformfactor(P)/sqrt(5)

Celestial body's `secondzonalharmonic` in spherical approximation (dimensionless).

```Julia
julia> secondzonalharmonic(Earth)
$(secondzonalharmonic(Earth))
```
"""
@pure secondzonalharmonic(P::Planet) = -dynamicformfactor(P)/sqrt(5)

@pure function _gravity(ϕ,P::Planet=Earth,U::US=Metric)
    β = latitudeparametric(ϕ,P)
    m,E = normal(oblateness(P)),lineareccentricity(P,U)
    a,b = semimajor(P,U),semiminor(P,U)
    q = m*eccentricity2(P)*q01(P)/3q0(P)
    sβ,cβ = sin(β),cos(β)
    g = gravitation(P,U)/(a*sqrt((a*sβ)^2+(b*cβ)^2))
    return g*((1+q)*sβ^2 + (1-m-q/2)*cβ^2)
end # geodetic gravity normal

"""
    gravity(ϕ,P::Planet)

Calculate total `gravity` at geodetic latitude `ϕ` on `Planet` surface (m⋅s⁻²).

```Julia
julia> gravity(0,Earth)
$(gravity(0,Earth))

julia> gravity(π/2,Earth)
$(gravity(π/2,Earth))
```
"""
@pure function gravity(ϕ::Real,P::Planet{f},U::US=Metric) where f
    sϕ2,ge,gp = sin(ϕ)^2,_gravity(0,P,U),_gravity(π/2,P,U)
    ge*((1+(normal(aspectratio(P))*(gp/ge)-1)*sϕ2)/sqrt(1-(f*(2-f))*sϕ2))
end # somigliana(1.0111032235724*π/4)

"""
    gravitygeodetic(h,ϕ,P::Planet)

Calculate total `gravity` at geodetic latitude `ϕ` and altitude `h` (m⋅s⁻²).
"""
@pure function gravitygeodetic(h,ϕ,P::Planet=Earth,U::US=Metric)
    ha,f,m = normal(h/semimajor(P,U)),flattening(P),oblateness(P)
    gravity(ϕ,P,U)*(1-2*(1+f+m-2f*sin(ϕ)^2)*ha+3*ha^2)
end

"""
    gravitycomponents(h,θ,P::Planet)

Calculate components of `gravity` at geocentric latitude `θ` and altitude `h` (m⋅s⁻²).
"""
@pure function gravitycomponents(h,θ,P::Planet=Earth,U::US=Metric)
    r = radius(θ,P,U)+h
    J2ar = 3*dynamicformfactor(P)*(semimajor(P,U)/r)^2
    sθ,cθ = sin(θ),cos(θ)
    g,ω = gravitation(P,U)/r^2,angularfrequency(P,U)
    Values(g*J2ar*sθ*cθ+r*ω^2*sθ*cθ,g*(1-J2ar/2*(3sθ^2-1))-r*(ω*cθ)^2)
end

@pure _gravity(h,θ,P::Planet=Earth,U::US=Metric) = norm(gravitycomponents(h,θ,P,U))

"""
    gravity(h,θ,P::Planet)

Calculate total `gravity` at geocentric latitude `θ` and altitude `h` (m⋅s⁻²).
"""
@pure function gravity(h,θ,P::Planet=Earth,U::US=Metric)
    gp,gp0 = _gravity(π/2,P,U),_gravity(0,π/2)
    _gravity(h,θ,P,U)*(1+((gp-gp0)/3gp)*sin(θ)^2)
end

if usingSimilitude
struct Atmosphere{n,P,U}
    a::Quantities{lapserate,U,n,Float64}
    h::Quantities{L,U,n,Float64} # altitude
    m::Values{n,Float64} # molar rate
    Atmosphere{P,U}(a::Values{n},h::Values{n},m::Values{n}) where {n,P,U} = new{n,P,normal(U)}(Quantities(lapserate,U,a),Quantities(L,U,h),m)
end
else
struct Atmosphere{n,P,U}
    a::Values{n,Float64} # lapse rate
    h::Values{n,Float64} # altitude
    m::Values{n,Float64} # molar rate
    Atmosphere{P,U}(a::Values{n},h::Values{n},m::Values{n}) where {n,P,U} = new{n,P,U}(a,h,m)
end
end

@doc """
    Atmosphere{n,P,U}

Temperature column of `P::Planet` with `U::UnitSystem` having `n` thermal layers.
""" Atmosphere

@pure Planet(::Atmosphere{n,P}) where {n,P} = P
@pure units(::Atmosphere{n,P,U}) where {n,P,U} = Quantity(U)
Atmosphere(a::Values{n},h::Values{n}) where n = Atmosphere{Earth}(a,h)
Atmosphere{P}(a::Values{n},h::Values{n}) where {n,P} = Atmosphere{P,Metric}(a,h)
Atmosphere{P,U}(a::Values{n},h::Values{n}) where {n,P,U} = Atmosphere{P,U}(a,h,zeros(Values{n}))
(U::UnitSystem)(A::Atmosphere{n,P,S}) where {n,P,S} = Atmosphere{P,U}(lapserate.(A.a,Ref(U),Ref(S)),length.(A.h,Ref(U),Ref(S)))

function display(A::Atmosphere)
    println(typeof(A))
    println(" a = ",A.a)
    println(" h = ",A.h)
    println(" m = ",A.m)
end

# local weather models

if usingSimilitude
struct Weather{ϕ,f,n,P,U} # ϕ ↦ radius, acceleration
    A::Atmosphere{n,P,U} # altitude lapse rate
    T::Quantities{Θ,U,n,Float64} # temperature
    p::Quantities{pressure,U,n,Float64} # pressure
    ρ::Quantities{density,U,n,Float64} # density
    Tc::Quantity{Θ,U,Float64}
    ha::Float64
    Weather{ϕ,f}(A::Atmosphere{n,P,U},T::Values{n},p::Values{n},ρ::Values{n},Tc,ha) where {ϕ,f,n,P,U} = new{ϕ,f,n,P,normal(U)}(A,Quantities(Θ,U,T),Quantities(pressure,U,p),Quantities(density,U,ρ),Quantity(Θ,U,Tc),ha)
end
(U::UnitSystem)(W::Weather{ϕ,f,n,P,S}) where {ϕ,f,n,P,S} = Weather{ϕ,f}(U(W.A),U(W.T),U(W.p),U(W.ρ))
else
struct Weather{ϕ,f,n,P,U} # ϕ ↦ radius, acceleration
    A::Atmosphere{n,P,U} # altitude lapse rate
    T::Values{n,Float64} # temperature
    p::Values{n,Float64} # pressure
    ρ::Values{n,Float64} # density
    Tc::Float64
    ha::Float64
    Weather{ϕ,f}(A::Atmosphere{n,P,U},T::Values{n},p::Values{n},ρ::Values{n},Tc,ha) where {ϕ,f,n,P,U} = new{ϕ,f,n,P,U}(A,T,p,ρ,Tc,ha)
end
(U::UnitSystem)(W::Weather{ϕ,f,n,P,S}) where {ϕ,f,n,P,S} = Weather{ϕ,f}(U(W.A),temperature.(W.T,Ref(U),Ref(S)),pressure.(W.p,Ref(U),Ref(S)),density.(W.ρ,Ref(U),Ref(S)))
end

@doc """
    Weather{ϕ,f,n,P,U}

Thermodynamic column state of fluid `f` at geodetic latitude `ϕ`, having `n` thermal `Atmosphere` layers (of `P::Planet` with `U::UnitSystem`).
Induces derived values `fluid`, `temperature`, `pressure`, `density`, `specificvolume`, `kinematic`, `heatcapacity`, `thermaldiffusivity`, `elasticity`, `specificimpedance`, `intensity`, `specificweight`, `geopotential`, and inherited values associated with `f::AbstractMole` derivations.
""" Weather

function Weather{ϕ}(A::Atmosphere{n,P,u},F::FluidState{f}) where {ϕ,f,n,P,u}
    T = zeros(Variables{n,Float64})
    p = zeros(Variables{n,Float64})
    ρ = zeros(Variables{n,Float64})
    μ = zeros(Variables{n,Float64})
    k = zeros(Variables{n,Float64})
    c = zeros(Variables{n,Float64})
    Δμ = zeros(Variables{n,Float64})
    Δk = zeros(Variables{n,Float64})
    Δc = zeros(Variables{n,Float64})
    U = Quantity(u)
    T0,p0 = normal(temperature(F,U)),normal(pressure(F,U))
    R,γ = normal(gasconstant(F,U)),normal(heatratio(F,U))
    r,g = normal(radiusgeodetic(ϕ,P,U)),normal(gravity(ϕ,P,U))
    T[1] = T0; p[1] = p0; ρ[1] = p0/(R*T0)
    μ[1],k[1],c[1] = normal(viscosity(F,U)),normal(thermalconductivity(F,U)),normal(heatvolume(F,U))
    Tc,ha = 0.0,0.0
    for i ∈ 2:n
        Δh = normal(A.h[i]-A.h[i-1])
        #T[i] = T[i-1]+A.a[i-1]*Δh
        T[i] = if isinf(normal(A.a[i-1]))
            Tc = (normal(A.a[i])*Δh*T[i]+T[i-1]^2-T[i]^2)/(normal(A.h[i])*Δh+2T[i-1]-2T[i])
            ha = Δh*(T[i-1]-Tc)/sqrt((T[i-1]-Tc)^2-(T[i]-Tc)^2)
            Tc+(T[i-1]-Tc)*sqrt(1-(Δh/ha)^2)
        else
            T[i-1]+normal(A.a[i-1])*Δh
        end
        vp,vρ = if iszero(normal(A.a[i-1]))
            val = exp((-g/R)*Δh/T[i])
            val,val
        else
            t,gRa = T[i]/T[i-1],-g/R/normal(A.a[i-1])
            t^gRa,t^(gRa-1)
        end
        p[i] = p[i-1]*vp
        ρ[i] = ρ[i-1]*vρ
        μ[i] = normal(viscosity(T[i],f,U))
        k[i] = normal(thermalconductivity(T[i],f,U))
        c[i] = normal(heatvolume(T[i],f,U))
        Δμ[i-1] = (μ[i]-μ[i-1])/Δh
        Δk[i-1] = (k[i]-k[i-1])/Δh
        Δc[i-1] = (c[i]-c[i-1])/Δh
        if i == n
            TT = T[end]+normal(A.a[end])*Δh
            Δμ[end] = (normal(viscosity(TT,f,U))-μ[i])/Δh
            Δk[end] = (normal(thermalconductivity(TT,f,U))-k[i])/Δh
            Δc[end] = (normal(heatvolume(TT,f,U))-c[i])/Δh
        end
    end
    Weather{ϕ,f}(A,Values(T),Values(p),Values(ρ),Tc,ha)
end

(A::Atmosphere)(F::FluidState=Air(288.16,atm,US(A)),r=6.356766e6,g=g₀) = Weather{ϕ}(A,F) # check this??
(A::Atmosphere)(T,p=atm,ϕ=1.0111032235724*π/4) = Weather{ϕ}(A,Air(T,p,units(A)))
(W::Weather)(h::Real=0) = (hG=normal(altgeopotent(h,W)); W(hG,layer(hG,W)))
function (W::Weather)(hG::Real,i)
    T = temperature(hG,i,W)
    FluidState{fluid(W),units(W)}(T,pressure(hG,T,i,W))
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
    T*temperature(S,U),a*lapserate(S,U),h*length(S,U),p*pressure(S,U),ρ*density(S,U)
end
@pure @inline Base.getindex(W::Weather,::Val{i}) where i = getindex(W,i)
@pure @inline lapserate(h::Real,W::Weather=Standard) = W.A.a[layer(h,W)]
@pure @inline layer(h::Real,W::Weather=Standard) = h≤normal(W.A.h[1]) ? 1 : (i=findfirst(x->normal(x)≥h,W.A.h); isnothing(i) ? length(W.A.h) : i-1)
@pure @inline layer(h::Real,W::Weather,U::US) = layer(length(h,units(W),U),W)

@pure Planet(::Weather{ϕ,f,n,P}) where {ϕ,f,n,P} = P
@pure units(::Weather{ϕ,f,n,P,U}) where {ϕ,f,n,P,U} = Quantity(U)
@pure fluid(::Weather{ϕ,f}=Standard) where {ϕ,f} = f

"""
    latitude(::Weather)

Geodetic latitude `ϕ` at `Weather` column location (rad).
"""
@pure latitude(W::Weather{ϕ}=Standard) where ϕ = ϕ

"""
    radius(::Weather)

Sea level radius `r` to planet's geodetic focus at `Weather` column location (m or ft).
"""
@pure radius(W::Weather=Standard,U::US=US(W)) = radiusgeodetic(latitude(W),Planet(W),U)

"""
    gravity(::Weather)

Sea level gravitational acceleration `g` at `Weather` column location (m⋅s⁻² or ft⋅s⁻²).
"""
@pure gravity(W::Weather=Standard,U::US=units(W)) = gravity(latitude(W),Planet(W),U)

@pure molecularmass(W::Weather=Standard,U::US=units(W)) = molecularmass(fluid(W),U)
@pure gasconstant(W::Weather=Standard,U::US=units(W)) = gasconstant(fluid(W),U)

# hG = geopotential altitude, h = geometric altitude

"""
    altabs(h::Real,W::Weather=Standard) = radius(W)+h

Absolute altitude from planet's center of gravity (m or ft).
"""
@pure @inline altabs(h::Real=0,W::Weather=Standard,U::US=units(W)) = radius(W,U)+Quantity(L,U,h)
@pure @inline altabs(h::Real,W::Weather,U::US,S::US) = altabs(h*length(S,U),W,U)

"""
    altgeopotent(h::Real,W::Weather=Standard) = h*radius(W)/altabs(h,W)

Geopotential altitude `hG` conversion from geometric altitude (m or ft).
"""
@pure altgeopotent(h::Real,W::Weather=Standard,U::US=US(W)) = (h/altabs(h,W,U))radius(W,U)
@pure altgeopotent(h::Real,W::Weather,U::US,S::US) = altgeopotent(h*length(S,U),W,U)

"""
    altgeometric(hG::Real,W::Weather=Standard) = radius(W)/(radius(W)/hG-1)

Geometric altitude `h` conversion from geopotential altitude (m or ft).
"""
@pure altgeometric(hG::Real,W::Weather=Standard,U::US=US(W)) = (r=radius(W,U); r/(r/hG-1))
@pure altgeometric(hG::Real,W::Weather,U::US,S::US) = altgeometric(hG*length(S,U),W,U)

"""
    gravity(h::Real=0,W::Weather=Standard) = gravity(W)*radius(W)^2/altabs(h,W)^2

Gravitational acceleration `g` at altitude `h` of `Weather` column (m⋅s⁻² or ft⋅s⁻²).
"""
@pure function gravity(h::Real,W::Weather=Standard,U::US=US(W))
    if h ≤ normal(0.007radius(W))
        (gravity(W,U)*radius(W,U)^2)/altabs(h,W,U)^2
    else
        gravitygeodetic(h,latitude(W),Planet(W),U)
    end
end
@pure gravity(h::Real,W::Weather,U::US,S::US) = gravity(h*length(S,U),W,U)

"""
    temperature(h::Real=0,::Weather=Standard)

Absolute temperature `T` at geometric altitude `h` of `Weather` location (K or °R).
"""
@pure function temperature(hG::Real,i,W::Weather=Standard,U::US=units(W))
    T0,a0,h0 = normal.(W[i,U])
    Quantity(Θ,U,if isinf(a0) # 1976 upper atmosphere
        Δh = hG-h0
        if a0 < 0 # 1976 elliptic layer
            W.Tc+(T0-W.Tc)*sqrt(1-(Δh/W.ha)^2)
        else # 1976 exponential layer
            r = radius(Earth1976)
            ξ = Δh*((r+h0)/(r+hG))
            1000-(1000-T0)*exp((-0.012/(1000-T0))*ξ)
        end
    else # standard
        iszero(a0) ? T0 : T0+a0*(hG-h0)
    end)
end

# k, μ, cᵥ, cₚ, γ, Pr, a, e, h

for op ∈ Intrinsic
    @eval @pure $op(hG::Real,i,W::Weather=Standard,U::US=US(W)) = $op(normal(temperature(hG,i,W,U)),fluid(W),U)
    #@eval @pure $op(h::Real,i,W::Weather,U::US,S::US) = $op(length(h,U,S),i,W,U)
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
    onreedom(h::Real=0,W::Weather) = freedom(temperature(h,W),fluid(W))

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
@pure pressure(hG::Real,i,W::Weather=Standard,U::US=units(W)) = pressure(hG,temperature(hG,i,W,U),i,W,U)
@pure function pressure(hG::Real,T::Number,i,W::Weather=Standard,U::US=units(W))
    g,R = gravity(W,U),gasconstant(W,U)
    T0,a,h0,p = W[i,U]
    p*normal(if iszero(a)
        exp((-g/R)*(Quantity(L,U,hG)-h0)/T)
    else
        (T/T0)^normal((-g/R)/a)
    end)
end

"""
    density(h::Real=0,W::Weather=Standard) = density(W(h))

Inertial mass per volume `ρ` at altitude `h` of `Weather` location (kg⋅m⁻³ or slugs⋅ft⁻³).
"""
@pure density(hG::Real,i,W::Weather=Standard,U::US=units(W)) = density(hG,temperature(hG,i,W,U),i,W,U)
#@pure density(hG::Real,T,i,W::Weather,U::US,S::US) = density(length(h,U,S),temperature(T,U,S),i,W,U)
@pure function density(hG::Real,T::Number,i,W::Weather=Standard,U::US=units(W))
    g,R = gravity(W,U),gasconstant(W,U)
    T0,a,h0,_,ρ = W[i,U]
    ρ*normal(if iszero(a)
        exp((-g/R)*(Quantity(L,U,hG)-h0)/T)
    else
        (T/T0)^normal((-g/R)/a-1)
    end)
end

"""
    kinematic(h::Real=0,W::Weather=Standard) = kinematic(W(h))

Kinematic viscosity ratio `ν` at altitude `h` of `Weather` location (m²⋅s⁻¹ or ft²⋅s⁻¹).
"""
@pure kinematic(hG::Real,i,W::Weather=Standard,U::US=units(W)) = kinematic(hG,temperature(hG,i,W,U),i,W,U)
#@pure kinematic(hG::Real,T,i,W::Weather,U::US,S::US) = kinematic(length(h,U,S),temperature(T,U,S),i,W,U)
@pure kinematic(hG::Real,T,i,W::Weather=Standard,U::US=units(W)) = viscosity(normal(T),fluid(W),U)/density(hG,T,i,W,U)

"""
    heatcapacity(h::Real=0,W::Weather=Standard) = heatcapacity(W(h))

Specific heat per mass at altitude `h` of `Weather` location (J⋅m⁻³⋅K⁻¹ or lb⋅ft⁻²⋅°R⁻¹).
"""
@pure function heatcapacity(hG::Real,i,W::Weather=Standard,U::US=units(W))
    T = normal(temperature(hG,i,W,U))
    heatpressure(T,fluid(W),U)*density(hG,T,i,W,U)
end

"""
    thermaldiffusivity(h::Real=0,W::Weather=Standard) = thermaldiffusivity(W(h))

Thermal diffusivity `α` at altitude `h` of `Weather` location (m²⋅s⁻¹ or ft²⋅s⁻¹).
"""
@pure function thermaldiffusivity(hG::Real,i,W::Weather=Standard,U::US=units(W))
    F = fluid(W)
    T = normal(temperature(hG,i,W,U))
    thermalconductivity(T,F,U)/heatpressure(T,F,U)/density(hG,T,i,W,U)
end

"""
    elasticity(h::Real=0,W::Weather=Standard) = elasticity(W(h))

Bulk modulus of elasticity `B` at altitude `h` of `Weather` location (Pa or slug⋅ft⁻¹⋅s⁻²).
"""
@pure function elasticity(hG,i,W::Weather=Standard,U::US=units(W))
    T = normal(temperature(hG,i,W,U))
    heatratio(T,fluid(W),U)*pressure(hG,T,i,W,U)
end

"""
    specificimpedance(h::Real=0,W::Weather=Standard) = impedance(W(h))

Specific acoustic resistance at altitude `h` of `Weather` (kg⋅m⁻³⋅s⁻¹ or slug⋅ft⁻³⋅s⁻¹).
"""
@pure function specificimpedance(hG::Real,i,W::Weather=Standard,U::US=units(W))
    T = normal(temperature(hG,i,W,U))
    density(hG,T,i,W,U)*sonicspeed(T,fluid(W),U)
end

"""
    intensity(h::Real=0,W::Weather=Standard) = intensity(W(h))

Instantaneous intensity `I` at altitude `h` of `Weather` at location (W⋅m⁻² or slug⋅s⁻³).
"""
@pure function intensity(hG::Real,i,W::Weather=Standard,U::US=units(W))
    T = temperature(hG,i,W,U)
    g,R = gravity(W,U),gasconstant(W,U)
    T0,a,h0,p,ρ = W[i,U]
    (p^2/ρ)*normal(if iszero(a)
        exp((-g/R)*(Quantity(L,U,hG)-h0)/T)
    else
        t,gRa = T/T0,normal((-g/R)/a)
        t^2gRa/t^(gRa-1)
    end)/sonicspeed(T,fluid(W),U)
end

# Grashof number

#=@pure function grashof(hG::Real,i,W::Weather=Standard,U::US=units(W))
    T = _temperature(hG,i,W,U)
    gravity(h,W,U)*(temperature(W,U)-T)*(h^3)/(T*kinematic(hG,T,i,W,U)^2)
end=#

"""
    specificweight(h::Real=0,W::Weather=Standard) = density(h,W)*gravity(h,W)

Specific weight at altitude `h` of `Weather` location (kg⋅m⁻²⋅s⁻² or slugs⋅ft⁻²⋅s⁻²).
"""
@pure specificweight(hG::Real,i,W::Weather=Standard,U::US=US(W)) = density(hG,i,W,U)*gravity(hG,W,U)

"""
    specificvolume(h::Real=0,W::Weather=Standard) = specificvolume(W(h))

Specific volume per mass `v` at altitude `h` of `Weather` location (m³⋅kg⁻¹, ft³⋅slug⁻¹).
"""
@pure specificvolume(hG::Real,i,W::Weather=Standard,U::US=US(W)) = inv(density(hG,i,W,U))

"""
    geopotential(h::Real=0,W::Weather=Standard) = gravity(h,W)*h

Specifc gravitational potential energy `g` at altitude `h` of `Weather` (m²⋅s⁻², ft²⋅s⁻²).
"""
@pure geopotential(h::Real,W::Weather=Standard,U::US=US(W)) = gravity(h,W,U)*h
@pure geopotential(h::Real,W::Weather,U::US,S::US) = geopotential(h*length(S,U),W,U)
@pure geopotential(W::Weather=Standard,U::US=Metric) = geopotential(0,W,U)

# common interface

for op ∈ (:temperature,:pressure,:density,:specificweight,:specificvolume,:specificimpedance,:thermaldiffusivity,:intensity,:heatcapacity,:kinematic,:elasticity,Intrinsic...) # grashof
    opratio = Symbol(op,:ratio)
    @eval begin
        export $op
        @pure function $op(h::Real,W::Weather=Standard,U::US=US(W))
            hG = normal(altgeopotent(h,W,U))
            $op(hG,layer(hG,W,U),W,U)
        end
        @pure $op(W::Weather=Standard,U::US=Metric) = $op(0,W,U)
        @pure $op(h::Real,W::Weather,U::US,S::US) = $op(length(h,U,S),W,U)
        #@pure $op(h,i,W::Weather,U::US,S::US) = $op(length(h,U,S),i,W,U)
    end
    op ∉ (:heatcapacity,:grashof,:geopotential) && @eval begin
        export $opratio
        @pure function $opratio(h::Real,W::Weather=Standard,U::US=US(W))
            hG = normal(altgeopotent(h,W,U))
            $opratio(hG,layer(hG,W,U),W,U)
        end
        @pure $opratio(hG::Real,i,W::Weather=Standard,U::US=US(W)) = $op(hG,i,W,U)/$op(W,U)
        @pure $opratio(h::Real,W::Weather,U::US,S::US) = $opratio(length(h,U,S),W,U)
        #@pure $opratio(h::Real,i,W::Weather,U::US,S::US) = $opratio(length(h,U,S),i,W,U)
    end
end

# data

include("planets.jl")

end # module
