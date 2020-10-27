module Geophysics

#   This file is part of Geophysics.jl. It is licensed under the MIT license
#   Geophysics Copyright (C) 2020 Michael Reed

using AbstractTensors
import Base: @pure

export Fluid, Atmosphere, Weather, layer, heatratio, gasconstant
export gravity, temperature, pressure, density, sonicspeed

# fluid models

struct Fluid{R,γ} end

@pure heatratio(::Fluid{R,γ}) where {R,γ} = γ
@pure gasconstant(::Fluid{R}) where R = R

struct Atmosphere{f,n}
    a::Values{n,Float64} # lapse rate
    h::Values{n,Float64} # altitude
    Atmosphere{f}(a::Values{n},h::Values{n}) where {f,n} = new{f,n}(a,h)
end

# local weather models

struct Weather{r,g,f,n} # radius, acceleration
    A::Atmosphere{f,n} # altitude lapse rate
    T::Values{n,Float64} # temperature
    p::Values{n,Float64} # pressure
    ρ::Values{n,Float64} # density
    Weather{r,g}(A::Atmosphere{f},T::Values{n},p::Values{n},ρ::Values{n}) where {r,g,f,n} = new{r,g,f,n}(A,T,p,ρ)
end

function Weather{r,g}(T0,p0,A::Atmosphere{f,n}) where {r,g,f,n} T = zeros(Variables{n,Float64})
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
    Weather{r,g}(A,Values(T),Values(p),Values(ρ))
end

(A::Atmosphere)(r,g,T,p) = Weather{r,g}(T,p,A)

@pure @inline Base.getindex(W::Weather,i::Int) = W.T[i],W.A.a[i],W.A.h[i],W.p[i],W.ρ[i]
@pure @inline Base.getindex(W::Weather,::Val{i}) where i = getindex(W,i)
@pure @inline lapserate(h::Real=0,W::Weather=Standard) = W[layer(h,W)]
@pure layer(h::Real=0,W::Weather=Standard) = h≤W.A.h[1] ? 1 : (i=findfirst(x->x≥h,W.A.h); isnothing(i) ? length(W.A.h) : i-1)

@pure fluid(::Weather{r,g,f}=Standard) where {r,g,f} = f
@pure heatratio(W::Weather=Standard) = heatratio(fluid(W))
@pure gasconstant(W::Weather=Standard) = gasconstant(fluid(W))
@pure radius(::Weather{r}=Standard) where r = r
@pure gravity(::Weather{r,g}=Standard) where {r,g} = g

# h = geopotential altitude, hG = geometric altitude

@pure @inline altabs(hG::Real=0,W::Weather=Standard) = radius(W)+hG
@pure @inline altratio(hG::Real=0,W::Weather=Standard) = radius(W)/altabs(hG,W)
@pure altgeopotent(hG::Real=0,W::Weather=Standard) = hG*altratio(hG,W)
@pure altgeometric(h::Real=0,W::Weather=Standard) = (r=radius(W); r/(r/h-1))
@pure gravity(hG::Real,W::Weather=Standard) = (gravity(W)*radius(W)^2)/altabs(hG,W)^2

# T = temperature

@pure temperature(h::Real=0,W::Weather=Standard) = temperature(h,layer(h,W),W)
@pure temperature(h::Real,i,W::Weather=Standard) = _temperature(altgeopotent(h,W),i,W)
@pure function _temperature(h::Real,i,W::Weather=Standard)
    T0,a0,h0 = W[i]
    return T0+a0*(h-h0)
end

# p = pressure

pressure(hG::Real=0,W::Weather=Standard) = pressure(hG,layer(hG,W),W)
function pressure(hG::Real,i,W::Weather=Standard)
    g,R = gravity(W),gasconstant(W)
    T0,a,h0,p = W[i]; h = altgeopotent(hG,W); T = _temperature(h,i,W)
    if iszero(a)
        p*exp((-g/R)*(h-h0)/T)
    else
        p*(T/T0)^((-g/R)/a)
    end
end

pressureratio(hG::Real=0,W::Weather=Standard) = pressureratio(hG,layer(hG,W),W)
pressureratio(hG::Real,i,W::Weather=Standard) = pressure(hG,i,W)/pressure(W)

# ρ = density

density(hG::Real=0,W::Weather=Standard) = density(hG,layer(hG,W),W)
function density(hG::Real,i,W::Weather=Standard)
    g,R = gravity(W),gasconstant(W)
    T0,a,h0,_,ρ = W[i]; h = altgeopotent(hG,W); T = _temperature(h,i,W)
    if iszero(a)
        ρ*exp((-g/R)*(h-h0)/T)
    else
        ρ*(T/T0)^((-g/R)/a-1)
    end
end

densityratio(hG::Real=0,W::Weather=Standard) = densityratio(hG,layer(hG,W),W)
densityratio(hG::Real,i,W::Weather=Standard) = density(hG,i,W)/density(W)

# speed of sound

@pure sonicspeed(hG::Real=0,W::Weather=Standard) = sonicspeed(hG,layer(hG,W),W)
@pure sonicspeed(hG::Real,i,W::Weather=Standard) = sqrt(heatratio(W)*gasconstant(W))*sqrt(temperature(hG,i,W))

include("planets.jl")

end # module
