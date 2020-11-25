#   This file is part of Geophysics.jl. It is licensed under the AGPL license
#   Geophysics Copyright (C) 2020 Michael Reed

export slug, ft, KJ1990, KJ2014, RK1990, RK2014
export mass, slugs, kilograms, poundal, meters, feet, rankine, kelvin, moles, molecules
export UnitSystem, CGS, CGS2019, Metric, SI2019, CODATA, Conventional, English
export Planck, PlanckGauss, Stoney, Hartree, Rydberg, Schrodinger, Electronic, Natural, NaturalGauss, QCD, QCDGauss, QCDoriginal

# unit systems

"""
    UnitSystem{kB,ħ,𝘤,μ₀,mₑ}

Standardized for engineering based on fundamental constants: `kB` Boltzmann's constant, `ħ` reduced Planck's constant, `𝘤` speed of light, `μ₀` vacuum permeability, and `mₑ` electron rest mass.
Primarily the `Metric` SI unit system is used in addition to the historic `English` engineering unit system.
These constants induce derived values for `avogadro`, `boltzmann`, `universal`, `planck`, `planckreduced`, `lightspeed`, `planckmass`, `atomicmass`, `protonmass`, `electronmass`, `newton`, `einstein`, `permeability`, `permittivity`, `coulomb`, and
additional constants `stefan`, `radiationintensity`, `impedance`, `charge`, `magneton`, `conductance`, `faraday`, `magneticflux`, `josephson`, `klitzing`, `hartree`, `rydberg`, `bohr`, `bohrreduced`, and `molarmass`.

Additional reference `UnitSystem` variants `CGS`, `CGS2019`, `SI2019`, `CODATA`, `Conventional`; along with several natural atomic units based on the fine structure constant `1/αinv` and the gravitational coupling constant `αG` (`Planck`, `PlanckGauss`, `Stoney`, `Hartree`, `Rydberg`, `Schrodinger`, `Electronic`, `Natural`, `NaturalGauss`, `QCD`, `QCDGauss`, and `QCDoriginal`).
""" #`Rᵤ,mᵤ,σ,ħ,μ₀,ε₀,kₑ,𝘦,𝔉,RK,Z₀,G₀`
struct UnitSystem{kB,ħ,𝘤,μ,mₑ} end
@pure boltzmann(::UnitSystem{k}) where k = k
@pure planckreduced(::UnitSystem{k,h}) where {k,h} = h
@pure lightspeed(::UnitSystem{k,h,c}) where {k,h,c} = c
@pure permeability(::UnitSystem{k,h,c,μ}) where {k,h,c,μ} = μ
@pure electronmass(::UnitSystem{k,h,c,μ,m}) where {k,h,c,μ,m} = m
# ΔνCs:s⁻¹, c:m⋅s⁻¹, h:kg⋅m²⋅s⁻¹, kB:kg⋅m²⋅s⁻²⋅K⁻¹, NA:mol⁻¹, Kcd: cd⋅sr⋅s³⋅kg⁻¹⋅m⁻²

const atm,𝘤,lbm = 101325.0,299792458.,32.17404856 # lb-f to pdl
const slug,ft,rankine,kelvin = 0.45359237lbm,𝘤/983571056.0,9/5,5/9
const kcalₜₕ,kcal₄,kcal₁₀,kcal₂₀,kcalₘ,kcalᵢₜ = 4184,4204,4185.5,4182,4190,4186.8
const calₜₕ,cal₄,cal₁₀,cal₂₀,calₘ,calᵢₜ = (kcalₜₕ,kcal₄,kcal₁₀,kcal₂₀,kcalₘ,kcalᵢₜ)./1e3
const kcal = kcalₜₕ; const cal = kcal/1000 # calₜₕ thermal calorie

@pure mass(U::UnitSystem,S::UnitSystem=Metric) = electronmass(U)/electronmass(S)
@pure mass(m::Real,U::UnitSystem=English,S::UnitSystem=Metric) = m*mass(U,S)

Base.display(U::UnitSystem) = println("UnitSystem{kB=$(boltzmann(U)),ħ=$(planckreduced(U)),𝘤=$(lightspeed(U)),μ₀=$(permeability(U)),mᵤ=$(electronmass(U))}")

# fundamental constants, αinv = (34259-1/4366.8123)/250 # 137.036 exactly?

const ΔνCs,Kcd,mP = 9192631770.0,683.0,2.176434e-8 # planck mass (kg)
const NA,kB,𝘩,𝘦 = 6.02214076e23,1.380649e-23,6.62607015e-34,1.602176634e-19
const μₑᵤ,μₚᵤ,αinv,R∞ = 1/1822.888486209,1.007276466621,137.035999084,10973731.5681601
const μ₀ = 2𝘩/𝘤/αinv/𝘦^2 # ≈ 4π*(1e-7+5.5e-17), exact charge
const ħ,δμ₀,μₚₑ,Rᵤ,mₑ = 𝘩/2π,μ₀-4π*1e-7,μₚᵤ/μₑᵤ,NA*kB,αinv^2*R∞*2𝘩/𝘤 # electron mass
const KJ1990,KJ2014,RK1990,RK2014 = 25812.807,25812.8074555,4.835979e14,4.835978525e14

# engineering units

const CGS = UnitSystem{1e10*Rᵤ*mₑ/μₑᵤ,1e7*ħ,100𝘤,4π,1000mₑ}()
const CGS2019 = UnitSystem{1e7*kB,1e7*ħ,100𝘤,1e7*μ₀,1000mₑ}()
const Metric = UnitSystem{Rᵤ*mₑ/μₑᵤ/0.001,ħ,𝘤,4π*1e-7,mₑ}()
const SI1976 = UnitSystem{8.31432mₑ/μₑᵤ/0.001,ħ,𝘤,4π*1e-7,mₑ}()
const SI2019 = UnitSystem{kB,ħ,𝘤,μ₀,mₑ}()
const CODATA = UnitSystem{Rᵤ*mₑ/μₑᵤ/0.001,2/KJ2014/RK2014^2/π,𝘤,2KJ2014/𝘤/αinv,mₑ}()
const Conventional = UnitSystem{1000Rᵤ*mₑ/μₑᵤ,2/KJ1990/RK1990^2/π,𝘤,2KJ1990/𝘤/αinv,mₑ}()
const English = UnitSystem{5.657302466e-24,ħ/slug/ft^2,𝘤/ft,4π,mₑ/slug}()

# natural units

const αG = (mₑ/mP)^2
const Planck = UnitSystem{1,1,1,1,√(4π*αG)}()
const PlanckGauss = UnitSystem{1,1,1,4π,√αG}()
const Stoney = UnitSystem{1,αinv,1,4π,√(αG*αinv)}()
const Hartree = UnitSystem{1,1,αinv,4π/αinv^2,1}()
const Rydberg = UnitSystem{1,1,2αinv,π/αinv^2,1/2}()
const Schrodinger = UnitSystem{1,1,αinv,4π/αinv^2,√αinv*mₑ/mP}()
const Electronic = UnitSystem{1,αinv,1,4π,1}()
const Natural = UnitSystem{1,1,1,1,1}()
const NaturalGauss = UnitSystem{1,1,1,4π,1}()
const QCD = UnitSystem{1,1,1,1,1/μₚₑ}()
const QCDGauss = UnitSystem{1,1,1,4π,1/μₚₑ}()
const QCDoriginal = UnitSystem{1,1,1,4π/αinv,1/μₚₑ}()

@pure molarmass(U::UnitSystem{1}) = 1
@pure molarmass(U::UnitSystem{boltzmann(CGS)}) = molarmass(Natural)
@pure molarmass(U::UnitSystem{kB}) where kB = molarmass(CGS)/1000
@pure molarmass(U::UnitSystem{1e7*kB}) = 1000molarmass(SI2019)
@pure molarmass(U::UnitSystem{kB}) = electronmass(U)*NA/μₑᵤ
@doc """
    molarmass(U) = avogadro(U)*electronmass(U)/μₑᵤ # 1/μₑᵤ = $(1/μₑᵤ)

Molar mass constant `Mᵤ` is ratio of the `molarmass` and `relativemass` of a chemical.
```Julia
julia> molarmass(CGS) # g⋅mol⁻¹
$(molarmass(CGS))

julia> molarmass(CGS2019) # g⋅mol⁻¹
$(molarmass(CGS2019))

julia> molarmass(Metric) # kg⋅mol⁻¹
$(molarmass(Metric))

julia> molarmass(SI2019) # kg⋅mol⁻¹
$(molarmass(SI2019))
```
""" molarmass, Mᵤ

@pure avogadro(U::UnitSystem) = μₑᵤ*molarmass(U)/electronmass(U)
@doc """
    avogadro(x) = universal(x)/boltzmann(x) # Mᵤ/atomicmass(x), Mᵤ ≈ 0.001-3.5e-13

Avogadro `NA` is `molarmass(x)/atomicmass(x)` number of atoms in a 12 g sample of C₁₂.
```Julia
julia> avogadro(SI2019) # mol⁻¹
$(avogadro(SI2019))

julia> avogadro(Metric) # mol⁻¹
$(avogadro(Metric))

julia> avogadro(English) # slug-mol⁻¹
$(avogadro(English))
```
""" avogadro, NA

# constants

@doc """
    planckreduced(x) = planck(x)/2π

Reduced Planck constant `ħ` is a Planck per radian (J⋅s⋅rad⁻¹ or ft⋅lb⋅s⋅rad⁻¹).

```Julia
julia> planckreduced(SI2019) # J⋅s⋅rad⁻¹
$(planckreduced(SI2019))

julia> planckreduced(SI2019)*lightspeed(SI2019) # J⋅m⋅rad⁻¹
$(planckreduced(SI2019)*lightspeed(SI2019))

julia> planckreduced(CODATA) # J⋅s⋅rad⁻¹
$(planckreduced(CODATA))

julia> planckreduced(Conventional) # J⋅s⋅rad⁻¹
$(planckreduced(Conventional))

julia> planckreduced(English) # ft⋅lb⋅s⋅rad⁻¹
$(planckreduced(English))
```
""" planckreduced, ħ

@pure planck(U::UnitSystem) = 2π*planckreduced(U)
@doc """
    planck(x) = 2π*planckreduced(x)

Planck constant `𝘩` is energy per electromagnetic frequency (J⋅s or ft⋅lb⋅s).

```Julia
julia> planck(SI2019) # J⋅s
$(planck(SI2019))

julia> planck(SI2019)*lightspeed(SI2019) # J⋅m
$(planck(SI2019)*lightspeed(SI2019))

julia> planck(CODATA) # J⋅s
$(planck(CODATA))

julia> planck(Conventional) # J⋅s
$(planck(Conventional))

julia> planck(English) # ft⋅lb⋅s
$(planck(English))
```
""" planck, 𝘩, hh

@doc """
    boltzmann(x) = universal(x)/avogadro(x)

Boltzmann constant `kB` is the entropy amount of a unit number microstate permutation.
```Julia
pressure*molecularmass == density*boltzmann*temperature
```
It satisfies the ideal gas law.

```Julia
julia> boltzmann(SI2019) # J⋅K⁻¹
$(boltzmann(SI2019))

julia> boltzmann(Metric) # J⋅K⁻¹
$(boltzmann(Metric))

julia> boltzmann(CGS) # erg⋅K⁻¹
$(boltzmann(CGS))

julia> boltzmann(SI2019)/planck(SI2019) # Hz⋅K⁻¹
$(boltzmann(SI2019)/planck(SI2019))

julia> boltzmann(SI2019)/calᵢₜ # calᵢₜ⋅K⁻¹
$(boltzmann(SI2019)/calᵢₜ)

julia> boltzmann(SI2019)/rankine/calᵢₜ # calᵢₜ⋅°R⁻¹
$(boltzmann(SI2019)/rankine/calᵢₜ)

julia> boltzmann(English) # ft⋅lb⋅°R⁻¹
$(boltzmann(English))

julia> boltzmann(SI2019)/planck(SI2019)/lightspeed(SI2019) # m⁻¹⋅K⁻¹
$(boltzmann(SI2019)/planck(SI2019)/lightspeed(SI2019))

julia> avogadro(SI2019)*boltzmann(SI2019)/calᵢₜ # calᵢₜ⋅mol⁻¹⋅K⁻¹
$(avogadro(SI2019)*boltzmann(SI2019)/calᵢₜ)

julia> 10log10(boltzmann(SI2019)/1) # dB(W⋅K⁻¹⋅Hz⁻¹)
$(10log10(boltzmann(SI2019)))
```
""" boltzmann, kB

@doc """
    lightspeed(U::UnitSystem) = 1/sqrt(μ₀*ε₀)

Speed of light in a vacuum `𝘤` for massless particles (m⋅s⁻¹ or ft⋅s⁻¹).

```Julia
julia> lightspeed(Metric) # m⋅s⁻¹
$(lightspeed(Metric))

julia> lightspeed(English) # ft⋅s⁻¹
$(lightspeed(English))
```
""" lightspeed, 𝘤, cc

@doc """
    permeability(U::UnitSystem) = 4π # Gaussian unit system

Magnetic permeability in a classical vacuum defined as `μ₀` in SI units (H⋅m⁻¹, kg⋅m²⋅C⁻²).

```Julia
julia> permeability(CGS) # abhenry⋅cm⁻¹
$(permeability(CGS))

julia> permeability(Metric) # H⋅m⁻¹
$(permeability(Metric))

julia> permeability(Conventional) # H⋅m⁻¹
$(permeability(Conventional))

julia> permeability(CODATA) # H⋅m⁻¹
$(permeability(CODATA))

julia> permeability(SI2019) # H⋅m⁻¹
$(permeability(SI2019))
```
""" permeability, μ₀, m0

@doc """
    electronmass(U::UnitSystem) = protonmass(U)/$μₚₑ # αinv^2*R∞*2𝘩/𝘤

Electron rest mass unit `mₑ` of subatomic particle with `-𝘦` elementary charge  (kg or slugs).
```Julia
julia> electronmass(Metric) # kg
$(electronmass(Metric))

julia> electronmass(Metric)/atomicmass(Metric) # Da
$μₑᵤ

julia> electronmass(Metric)*lightspeed(Metric)^2 # J
$(electronmass(Metric)*lightspeed(Metric)^2)

julia> electronmass(English) # slugs
$(electronmass(English))
```
""" electronmass, mₑ, me

@pure atomicmass(U::UnitSystem) = electronmass(U)/μₑᵤ
@doc """
    atomicmass(U::UnitSystem) = Mᵤ/avogadro(U) # $(molarmass(SI2019)) ≈ 0.001-3.5e-13

Atomic mass unit `mᵤ` of 1/12 of the C₁₂ carbon-12 atom's mass  (kg or slugs).
```Julia
julia> atomicmass(Metric) # kg
$(atomicmass(Metric))

julia> atomicmass(Metric)/electronmass(Metric) # mₑ
$(atomicmass(Metric)/electronmass(Metric))

julia> atomicmass(Metric)*lightspeed(Metric)^2 # J
$(atomicmass(Metric)*lightspeed(Metric)^2)

julia> atomicmass(English) # slugs
$(atomicmass(English))
```
""" atomicmass, mᵤ, mu

@pure protonmass(U::UnitSystem) =  μₚᵤ*atomicmass(U)
@doc """
    protonmass(U::UnitSystem) = $(μₚᵤ)atomicmass(U)

Proton mass unit `mₚ` of subatomic particle with `+𝘦` elementary charge  (kg or slugs).
```Julia
julia> protonmass(Metric) # kg
$(protonmass(Metric))

julia> protonmass(Metric)/atomicmass(Metric) # mᵤ
$(protonmass(Metric)/atomicmass(Metric))

julia> protonmass(Metric)/electronmass(Metric) # mₑ
$(protonmass(Metric)/electronmass(Metric))
```
""" protonmass, mₚ, mp

@pure planckmass(U::UnitSystem) = mass(mP,U)
@doc """
    planckmass(U::UnitSystem) = mass($(planckmass(Metric)),U)

Planck mass factor `mP` from the gravitational coupling constant `αG` (kg or slugs).
```Julia
juila> planckmass(Metric) # kg
$(planckmass(Metric))

juila> planckmass(Metric)/atomicmass(Metric) # mᵤ
$(planckmass(Metric)/atomicmass(Metric))

juila> planckmass(Metric)/sqrt(8π) # kg
$(planckmass(Metric)/sqrt(8π))
```
""" planckmass, mP

@pure newton(U::UnitSystem) = lightspeed(U)*planckreduced(U)/planckmass(U)^2
@doc """
    newton(U::UnitSystem) = lightspeed(U)*planckreduced(U)/planckmass(U)^2

Universal gravitational constant `GG` of Newton's law (m³⋅kg⁻¹⋅s⁻² or ft³⋅slug⁻¹⋅s⁻²).
```Julia
juila> newton(Metric) # m³⋅kg⁻¹⋅s⁻²
$(newton(Metric))

julia> newton(English) # ft³⋅slug⁻¹⋅s⁻²
$(newton(English))
```
""" newton, GG

@pure einstein(U::UnitSystem) = 8π*newton(U)/lightspeed(U)^4
@doc """
    einstein(U::UnitSystem) = 8π*newton(U)/lightspeed(U)^4

Einstein's gravitational constant from the Einstein field equations (? or ?).
```Julia
julia> einstein(Metric) # ?
$(einstein(Metric))
```
""" einstein, κ

@pure universal(U::UnitSystem) = boltzmann(U)*avogadro(U)
@doc """
    universal(x) = boltzmann(x)*avogadro(x)

Universal gas constant `Rᵤ` is factored into specific `gasconstant(x)*molarmass(x)` values.
```Julia
pressure*molarmass == density*universal*temperature
```
It satisfies the ideal gas law.

```Julia
julia> universal(SI2019) # J⋅K⁻¹⋅mol⁻¹
$(universal(SI2019))

julia> universal(English)*lbm/2116.2 # atm⋅ft³⋅R⁻¹⋅slug-mol⁻¹
$(universal(English)*lbm/2116.2)

julia> universal(Metric)/cal # cal⋅K⁻¹⋅mol⁻¹
$(universal(Metric)/cal)

julia> universal(Metric)/pressure(Earth1959) # atm⋅m³⋅K⁻¹⋅mol⁻¹
$(universal(Metric)/atm)

julia> universal(CGS) # erg⋅K⁻¹⋅mol⁻¹
$(universal(CGS))

julia> universal(English) # ft⋅lb⋅°R⁻¹⋅slug-mol⁻¹
$(universal(English))
```
The 1976 United States Standard Atmosphere used R* = 8.31432 exactly.
""" universal, Rᵤ, Ru

@pure stefan(U::UnitSystem) = 2π^5*boltzmann(U)^4/(15planck(U)^3*lightspeed(U)^2)
@doc """
    stefan(U::UnitSystem) = 2π^5*boltzmann(U)^4/(15planck(U)^3*lightspeed(U)^2)

Stefan-Boltzmann proportionality `σ` of black body radiation (W⋅m⁻²⋅K⁻⁴ or ?⋅ft⁻²⋅°R⁻⁴).

```Julia
julia> stefan(Metric) # W⋅m⁻²⋅K⁻⁴
$(stefan(Metric))

julia> stefan(CGS) # erg⋅cm⁻²⋅s⁻¹⋅K⁻⁴
$(stefan(CGS))

julia> stefan(Metric)*24*60^2/(cal*100^2) # cal⋅cm⁻²⋅day⁻¹⋅K⁻⁴
$(stefan(Metric)*24*0.6^2/cal)

julia> stefan(English) # lb⋅s⁻¹⋅ft⁻³⋅°R⁻⁴
$(stefan(English))
```
""" stefan, σ, SB

"""
    radiationdensity(U::UnitSystem) = 4stefan(U)/lightspeed(U)

Raditation density constant of black body radiation (J⋅m⁻³⋅K⁻⁴ or lb⋅ft⁻²⋅°R⁻⁴).

```Julia
julia> radiationdensity(Metric) # J⋅m⁻³⋅K⁻⁴
$(radiationdensity(Metric))

julia> radiationdensity(CGS) # erg⋅cm⁻³⋅K⁻⁴
$(radiationdensity(CGS))

julia> radiationdensity(English) # lb⋅ft⁻²⋅°R⁻⁴
$(radiationdensity(English))
```
"""
@pure radiationdensity(U::UnitSystem) = 4stefan(U)/lightspeed(U)

@pure permittivity(U::UnitSystem) = inv(permeability(U)*lightspeed(U)^2)
@doc """
    permittivity(U::UnitSystem) = 1/permeability(U)/lightspeed(U)^2

Dielectric permittivity constant `ε₀` of a classical vacuum (C²⋅N⁻¹⋅m⁻²).

```Julia
julia> permittivity(Metric) # F⋅m⁻¹
$(permittivity(Metric))

julia> permittivity(Conventional) # F⋅m⁻¹
$(permittivity(Conventional))

julia> permittivity(CODATA) # F⋅m⁻¹
$(permittivity(CODATA))

julia> permittivity(SI2019) # F⋅m⁻¹
$(permittivity(SI2019))
```
""" permittivity, ε₀, e0

@pure coulomb(U::UnitSystem) = inv(4π*permittivity(U))
@doc """
    coulomb(U::UnitSystem) = 1/4π/permittivity(U)

Electrostatic proportionality constant `kₑ` for the Coulomb's law force (N⋅m²⋅C⁻²).

```Julia
julia> coulomb(Metric) # N⋅m²⋅C⁻²
$(coulomb(Metric))

julia> coulomb(Metric)/lightspeed(Metric)^2 # (N·s²⋅C⁻²)⋅𝘤²
$(coulomb(Metric)/lightspeed(Metric)^2)

julia> coulomb(Conventional)/lightspeed(Conventional)^2 # (N·s²⋅C⁻²)⋅𝘤²
$(coulomb(Conventional)/lightspeed(Conventional)^2)

julia> coulomb(CODATA)/lightspeed(CODATA)^2 # (N·s²⋅C⁻²)⋅𝘤²
$(coulomb(CODATA)/lightspeed(CODATA)^2)

julia> coulomb(SI2019)/lightspeed(SI2019)^2 # (N·s²⋅C⁻²)⋅𝘤²
$(coulomb(SI2019)/lightspeed(SI2019)^2)
```
""" coulomb, kₑ, ke

@pure impedance(U::UnitSystem) = permeability(U)*lightspeed(U)
@doc """
    impedance(U::UnitSystem) = permeability(U)*lightspeed(U)

Vacuum impedance of free space `Z₀` is magnitude ratio of electric to magnetic field (Ω).
```Julia
julia> impedance(Metric) # Ω
$(impedance(Metric))

julia> impedance(Conventional) # Ω
$(impedance(Conventional))

julia> impedance(CODATA) # Ω
$(impedance(CODATA))

julia> impedance(SI2019) # Ω
$(impedance(SI2019))

julia> 120π # 3e8*μ₀ # Ω
$(120π)
```
""" impedance, Z₀, Z0

@pure charge(U::UnitSystem) = sqrt(2planck(U)/impedance(U)/αinv) # fine structure
@doc """
    charge(U::UnitSystem) = sqrt(2𝘩/$(αinv)impedance(U)) # faraday(U)/avogadro(U)

Quantized elementary charge `𝘦` of a proton or electron `2/(klitzing(U)*josephson(U))` (C).
```Julia
julia> charge(SI2019) # C
$(charge(SI2019))

julia> charge(Metric) # C
$(charge(Metric))

julia> charge(CODATA) # C
$(charge(CODATA))

julia> charge(Conventional) # C
$(charge(Conventional))

julia> 10lightspeed(Metric)*charge(metric) # statC
$(10lightspeed(Metric)*charge(Metric))

julia> charge(Planck) # sqrt(4π/αinv)
$(charge(Planck))
```
""" charge, 𝘦, ee

@pure faraday(U::UnitSystem) = charge(U)*avogadro(U)
@doc """
    faraday(U::UnitSystem) = charge(U)*avogadro(U)

Electric charge per mole of electrons `𝔉` based on elementary charge (C⋅mol⁻¹).
```Julia
julia> faraday(SI2019) # C⋅mol⁻¹
$(faraday(SI2019))

julia> faraday(Metric) # C⋅mol⁻¹
$(faraday(Metric))

julia> faraday(CODATA) # C⋅mol⁻¹
$(faraday(CODATA))

julia> faraday(Conventional) # C⋅mol⁻¹
$(faraday(Conventional))

julia> faraday(Metric)/kcal # kcal⋅(V-g-e)⁻¹
$(faraday(Metric)/kcal)

julia> faraday(Metric)/3600 # A⋅h⋅mol⁻¹
$(faraday(Metric)/3600)
```
""" faraday, 𝔉, FF

@pure josephson(U::UnitSystem) = 2charge(U)/planck(U)
@doc """
    josephson(U::UnitSystem) = 2charge(U)/planck(U) # 1/magneticflux(U)

Josephson constant `KJ` relating potential difference to irradiation frequency (Hz⋅V⁻¹).
```Julia
julia> josephson(SI2019) # Hz⋅V⁻¹
$(josephson(SI2019))

julia> josephson(Metric) # Hz⋅V⁻¹
$(josephson(Metric))

julia> josephson(Conventional) # Hz⋅V⁻¹
$(josephson(Conventional))

julia> josephson(CODATA) # Hz⋅V⁻¹
$(josephson(CODATA))
```
""" josephson, KJ

@pure magneticflux(U::UnitSystem) = inv(josephson(U))
@doc """
    magneticflux(U::UnitSystem) = planck(U)/2charge(U)

Magnetic flux quantum `Φ₀` is `1/josephson(U)` (Wb).
```Julia
julia> magneticflux(SI2019) # Wb
$(magneticflux(SI2019))

julia> magneticflux(Metric) # Wb
$(magneticflux(Metric))

julia> magneticflux(Conventional) # Wb
$(magneticflux(Conventional))

julia> magneticflux(CODATA) # Wb
$(magneticflux(CODATA))
```
""" magneticflux, Φ₀

@pure klitzing(U::UnitSystem) = planck(U)/charge(U)^2
@doc """
    klitzing(U::UnitSystem) = 2/conductance(U)

Quantized Hall resistance `RK` (Ω).
```Julia
julia> klitzing(SI2019) # Ω
$(klitzing(SI2019))

julia> klitzing(Metric) # Ω
$(klitzing(Metric))

julia> klitzing(Conventional) # Ω
$(klitzing(Conventional))

julia> klitzing(CODATA) # Ω
$(klitzing(CODATA))
```
""" klitzing, RK

@pure conductance(U::UnitSystem) = 2charge(U)^2/planck(U)
@doc """
    conductance(U::UnitSystem) = 2charge(U)^2/𝘩 # 2/klitzing(U)

Conductance quantum `G₀` is a quantized unit of electrical conductance (S).
```Julia
julia> conductance(SI2019) # S
$(conductance(SI2019))

julia> conductance(Metric) # S
$(conductance(Metric))

julia> conductance(Conventional) # S
$(conductance(Conventional))

julia> conductance(CODATA) # S
$(conductance(CODATA))
```
""" conductance, G₀, G0

@pure hartree(U::UnitSystem) = electronmass(U)*(lightspeed(U)/αinv)^2
@doc """
    hartree(U::UnitSystem) = electronmass(U)*(lightspeed(U)/$αinv)^2 # mₑ*(𝘤/αinv)^2

Hartree electric potential energy `Eₕ` of the hydrogen atom at ground state is `2R∞*𝘩*𝘤` (J).
```Julia
julia> hartree(Metric) # J
$(hartree(Metric))

julia> hartree(CGS) # erg
$(hartree(CGS))

julia> hartree(Metric)*avogadro(Metric)/1000 # kJ⋅mol⁻¹
$(hartree(Metric)*avogadro(Metric)/1000)

julia> hartree(Metric)*avogadro(Metric)/kcal # kcal⋅mol⁻¹
$(hartree(Metric)*avogadro(Metric)/kcal)

julia> 2rydberg(Metric)/100 # Eₕ/𝘩/𝘤/100 cm⁻¹
$(hartree(Metric)/planck(Metric)/lightspeed(Metric)/100)

julia> hartree(Metric)/planck(Metric)/10^12 # THz
$(hartree(Metric)/planck(Metric))

julia> hartree(Metric)/boltzmann(Metric) # K
$(hartree(Metric)/boltzmann(Metric))
```
In a Gaussian unit system where `4π*ε₀ == 1` the Hartree energy is `𝘦^2/a₀`.
""" hartree, Eₕ, Eh

@pure rydberg(U::UnitSystem) = hartree(U)/2planck(U)/lightspeed(U)
@doc """
    rydberg(U::UnitSystem) = hartree(U)/2planck(U)/lightspeed(U) # Eₕ/2𝘩/𝘤

Rydberg constant `R∞` is lowest energy photon capable of ionizing atom at ground state (m⁻¹).
```Julia
julia> rydberg(Metric) # m⁻¹
$(rydberg(Metric))
```
The Rydberg constant for hydrogen `RH` is `R∞*mₚ/(mₑ+mₚ)` (m⁻¹).
```Julia
julia> rydberg(Metric)*protonmass(Metric)/(electronmass(Metric)+protonmass(Metric))
$(rydberg(Metric)*protonmass(Metric)/(electronmass(Metric)+protonmass(Metric)))
```
Rydberg unit of photon energy `Ry` is `𝘩*𝘤*R∞` or `Eₕ/2` (J).
```Julia
julia> hartree(Metric)/2
$(hartree(Metric)/2)
```
Rydberg photon frequency `𝘤*R∞` or `Eₕ/2𝘩` (Hz)
```Julia
julia> lightspeed(Metric)*rydberg(Metric)
$(lightspeed(Metric)*rydberg(Metric))
```
Rydberg wavelength `1/R∞` (m)
```Julia
julia> 1/rydberg(Metric)
$(1/rydberg(Metric))

julia> 1/rydberg(Metric)/2π # angular wavelength
$(1/rydberg(Metric)/2π)
```
Precision measurements of the Rydberg constants are within a relative standard uncertainty of under 2 parts in 10¹², and is chosen to constrain values of other physical constants.
""" rydberg, R∞, RH, Ry

@pure plancklength(U::UnitSystem) = sqrt(planckreduced(U)*newton(U)/lightspeed(U)^3)
@pure bohr(U::UnitSystem) = αinv*planckreduced(U)/electronmass(U)/lightspeed(U)
@doc """
    bohr(U) = $αinv*planckreduced(U)/electronmass(U)/lightspeed(U)

Bohr radius of the hydrogen atom in its ground state `a₀` (m).
```Julia
julia> bohr(Metric) # m
$(bohr(Metric))

julia> 12bohr(English) # in
$(12bohr(English))

julia> bohr(Metric)/plancklength(Metric) # ℓP
$(bohr(Metric)/plancklength(Metric))
```
""" bohr, a₀, a0

"""
    bohrreduced(U::UnitSystem) = electronmass(U)/bohr(U)/$μₚₑ

Reduced Bohr radius including the effect of reduced mass in hydrogen atom (m).
```Julia
julia> bohrreduced(Metric) # m
$(bohrreduced(Metric))

julia> bohrreduced(Metric) # a₀
$(bohrreduced(Metric)/bohr(Metric))
```
"""
@pure bohrreduced(U::UnitSystem) = bohr(U)*(1+1/μₚₑ)

@pure electronradius(U::UnitSystem) = planckreduced(U)/electronmass(U)/lightspeed(U)/αinv
@doc """
    electronradius(U) = planckreduced(U)/electronmass(U)/lightspeed(U)/$αinv

Classical electron radius or Lorentz radius or Thomson scattering length (m).
```Julia
julia> electronradius(Metric) # m
$(electronradius(Metric))

julia> electronradius(CODATA) # m
$(electronradius(CODATA))

julia> electronradius(Conventional) # m
$(electronradius(Conventional))
```
""" electronradius, rₑ, re

@pure magneton(U::UnitSystem) = charge(U)*planckreduced(U)/2electronmass(U)
"""
    magneton(U::UnitSystem) = charge(U)*planckreduced(U)/2electronmass(U)

Bohr magneton `μB` natural unit for expressing magnetic moment of electron (J⋅T⁻¹).
```Julia
julia> magneton(SI2019) # J⋅T⁻¹
$(magneton(SI2019))

julia> magneton(Metric) # J⋅T⁻¹
$(magneton(Metric))

julia> magneton(CODATA) # J⋅T⁻¹
$(magneton(CODATA))

julia> magneton(Conventional) # J⋅T⁻¹
$(magneton(Conventional))

julia> magneton(CGS2019) # erg⋅T⁻¹
$(magneton(CGS2019))

julia> magneton(Hartree) # 𝘤⋅ħ⋅mₑ⁻¹
$(magneton(Hartree))
```
""" magneton, μB

# == Metric is different
const κ = einstein(SI2019)
const GG = newton(SI2019)
const σ = stefan(SI2019) #
const μB = magneton(SI2019) #
const ε₀ = permittivity(SI2019) #
const kₑ = coulomb(SI2019) #
const mₚ = protonmass(SI2019)
const mᵤ = atomicmass(SI2019)
const Mᵤ = molarmass(SI2019)
const 𝔉 = faraday(SI2019) #
const Φ₀ = magneticflux(SI2019) #
const Z₀ = impedance(SI2019) #
const G₀ = conductance(SI2019) #
const Eₕ = hartree(SI2019)
const a₀ = bohr(SI2019)
const rₑ = electronradius(SI2019)
const ℓP = plancklength(SI2019)
const RK = klitzing(SI2019) #
const KJ = josephson(SI2019) #
const RH,Ry = R∞*mₚ/(mₑ+mₚ),𝘩*𝘤*R∞
const Mu,Ru,SB,hh,cc,m0,e0,ke,me,mp,mu,ee,FF,Z0,G0,Eh,a0,re = Mᵤ,Rᵤ,σ,𝘩,𝘤,μ₀,ε₀,kₑ,mₑ,mₚ,mᵤ,𝘦,𝔉,Z₀,G₀,Eₕ,a₀,rₑ
export κ, GG, NA, kB, Rᵤ, σ, 𝘩, ħ, 𝘤, μ₀, ε₀, kₑ, mₑ, mₚ, mᵤ, 𝘦, 𝔉, Φ₀, Z₀, G₀, Eₕ, R∞, a₀, rₑ, KJ, RK, Ru, SB, hh, cc, m0, e0, ke, me, mp, mu, ee, FF, Z0, G0, Eh, a0, re
export αG, αinv, μₚₑ, μₑᵤ, μₚᵤ, mpe, meu, mpu, mP, δμ₀, Mᵤ, Mu, RH, Ry, ΔνCs, Kcd
export cal, kcal, calₜₕ, kcalₜₕ, calᵢₜ, kcalᵢₜ, SI, SI1976, ℓP, plancklength
const mpe, mea, mpu, SI = μₚₑ, μₑᵤ, μₚᵤ, SI2019

export electronmass, protonmass, atomicmass, planckmass, stefan, radiationdensity, einstein, impedance, charge, faraday, josephson, klitzing, hartree, rydberg, bohr, bohrreduced, electronradius, conductance, magneticflux, magneton, molarmass

const Constants = (:newton,:avogadro,:boltzmann,:planck,:planckreduced,:lightspeed,:universal,:permeability,:permittivity,:coulomb)

const Properties = (:units,:molarmass,:molecularmass,:gasconstant,Constants...)

const Intrinsic = (:viscosity,:conductivity,:heatvolume,:heatpressure,:heatratio,:prandtl,:sonicspeed,:freedom,:energy,:enthalpy)

# unit systems

@doc """
    Metric::UnitSystem{Rᵤ*mₑ/μₑᵤ/0.001,ħ,𝘤,4π*1e-7,mₑ}

Systeme International d'Unites (the SI units) adopted as the preffered `UnitSystem`.

```Julia
julia> boltzmann(Metric) # J⋅K⁻¹
$(boltzmann(Metric))

julia> planckreduced(Metric) # J⋅s⋅rad⁻¹
$(planckreduced(Metric))

julia> lightspeed(Metric) # m⋅s⁻¹
$(lightspeed(Metric))

julia> permeability(Metric) # H⋅m⁻¹
$(permeability(Metric))

julia> electronmass(Metric) # kg
$(electronmass(Metric))
```
""" Metric

@doc """
    English::UnitSystem{5.657302466e-24,ħ/slug/ft^2,𝘤/ft,4π,mₑ/slug}

Engineering `UnitSystem` historically used by Britain and United States.

```Julia
julia> boltzmann(English) # ft⋅lb⋅°R⁻¹
$(boltzmann(English))

julia> planckreduced(English) # ft⋅lb⋅s⋅rad⁻¹
$(planckreduced(English))

julia> lightspeed(English) # ft⋅s⁻¹
$(lightspeed(English))

julia> permeability(English) # slug⋅ft²⋅?⁻²
$(permeability(English))

julia> electronmass(English) # kg
$(electronmass(English))
```
""" English

# engineering units

@doc """
    CGS::UnitSystem{1e10*Rᵤ*mₑ/μₑᵤ,1e7*ħ,100𝘤,4π,1000mₑ}

Centimetre-gram-second `UnitSystem` variant of `Metric` system based on factors of `1e2,1e3`.

```Julia
julia> boltzmann(CGS) # erg⋅K⁻¹
$(boltzmann(CGS))

julia> planckreduced(CGS) # erg⋅s⋅rad⁻¹
$(planckreduced(CGS))

julia> lightspeed(CGS) # cm⋅s⁻¹
$(lightspeed(Metric))

julia> permeability(CGS) # erg⋅A⁻²⋅cm⁻¹
$(permeability(CGS))

julia> electronmass(CGS) # g
$(electronmass(CGS))
```
""" CGS

@doc """
    CGS2019::UnitSystem{1e7*kB,1e7*ħ,100𝘤,1e7*μ₀,1000mₑ}

Centimetre-gram-second `UnitSystem` variant of the tuned `SI2019` unit specification.

```Julia
julia> boltzmann(CGS2019) # erg⋅K⁻¹
$(boltzmann(CGS2019))

julia> planckreduced(CGS2019) # erg⋅s⋅rad⁻¹
$(planckreduced(CGS2019))

julia> lightspeed(CGS2019) # cm⋅s⁻¹
$(lightspeed(CGS2019))

julia> permeability(CGS2019) # erg⋅A⁻²⋅cm⁻¹
$(permeability(CGS2019))

julia> electronmass(CGS2019 # g
$(electronmass(CGS2019))
```
""" CGS2019

@doc """
    CODATA::UnitSystem{Rᵤ*mₑ/μₑᵤ/0.001,2/KJ2014/RK2014^2/π,𝘤,2KJ2014/𝘤/αinv,mₑ}

Metric `UnitSystem` based on Committee on Data of the International Science Council.

```Julia
julia> boltzmann(CODATA) # J⋅K⁻¹
$(boltzmann(CODATA))

julia> planckreduced(CODATA) # J⋅s⋅rad⁻¹
$(planckreduced(CODATA))

julia> lightspeed(CODATA) # m⋅s⁻¹
$(lightspeed(CODATA))

julia> permeability(CODATA) # H⋅m⁻¹
$(permeability(CODATA))

julia> electronmass(CODATA) # kg
$(electronmass(CODATA))
```
""" CODATA

@doc """
    SI2019::UnitSystem{kB,ħ,𝘤,μ₀,mₑ}

Systeme International d'Unites (the SI units) with `μ₀` for a tuned `charge` exactly.

```Julia
julia> boltzmann(SI2019) # J⋅K⁻¹
$(boltzmann(SI2019))

julia> planckreduced(SI2019) # J⋅s⋅rad⁻¹
$(planckreduced(SI2019))

julia> lightspeed(SI2019) # m⋅s⁻¹
$(lightspeed(SI2019))

julia> permeability(SI2019) # H⋅m⁻¹
$(permeability(CODATA))

julia> electronmass(SI2019) # kg
$(electronmass(SI2019))
```
""" SI2019, SI

@doc """
    Conventional::UnitSystem{Rᵤ*mₑ/μₑᵤ/0.001,2/KJ1990/RK1990^2/π,𝘤,2KJ1990/𝘤/αinv,mₑ}

Conventional electronic `UnitSystem` with 1990 tuned `josephson` and `klitzing` constants.

```Julia
julia> boltzmann(Conventional) # J⋅K⁻¹
$(boltzmann(Conventional))

julia> planckreduced(Conventional) # J⋅s⋅rad⁻¹
$(planckreduced(Conventional))

julia> lightspeed(Conventional) # m⋅s⁻¹
$(lightspeed(Conventional))

julia> permeability(Conventional) # H⋅m⁻¹
$(permeability(Conventional))

julia> electronmass(Conventional) # kg
$(electronmass(Conventional))
```
""" Conventional

textunits(U,S) = """
```Julia
julia> boltzmann($S)
$(boltzmann(U))

julia> planckreduced($S)
$(planckreduced(U))

julia> lightspeed($S)
$(lightspeed(U))

julia> permeability($S)
$(permeability(U))

julia> electronmass($S)
$(electronmass(U))
```
"""

@doc """
    Planck::UnitSystem{1,1,1,1,√(4π*αG)}

Planck `UnitSystem` with the `electronmass` value `√(4π*αG)` using gravitational coupling.

$(textunits(Planck,:Planck))
""" Planck

@doc """
    PlanckGauss::UnitSystem{1,1,1,4π,√αG}

Planck (Gauss) `UnitSystem` with `permeability` of `4π` and `electronmass` coupling `√αG`.

$(textunits(PlanckGauss,:PlanckGauss))
""" PlanckGauss

@doc """
    Stoney::UnitSystem{1,αinv,1,4π,√(αG*αinv)}

Stoney `UnitSystem` with `permeability` of `4π` and `electronmass` coupling `√(αG*αinv)`.

$(textunits(Stoney,:Stoney))
""" Stoney

@doc """
    Hartree::UnitSystem{1,1,αinv,4π/αinv^2,1}

Hartree atomic `UnitSystem` with `lightspeed` of `αinv` and `permeability` of `4π/αinv^2`.

$(textunits(Hartree,:Hartree))
""" Hartree

@doc """
    Rydberg::UnitSystem{1,1,2αinv,π/αinv^2,1/2}

Rydberg `UnitSystem` with `lightspeed` of `2αinv` and `permeability` of `π/αinv^2`.

$(textunits(Rydberg,:Rydberg))
""" Rydberg

@doc """
    Schrodinger::UnitSystem{1,1,αinv,4π/αinv^2,√αinv*mₑ/mP}

Schrodinger `UnitSystem` with `permeability` of `4π/αinv^2` and `electronmass` of `√αinv*mₑ/mP`.

$(textunits(Schrodinger,:Schrodinger))
""" Schrodinger

const Electronic = UnitSystem{1,αinv,1,4π,1}()
@doc """
    Electronic::UnitSystem{1,αinv,1,4π,1}

Electronic `UnitSystem` with `planckreduced` of `αinv` and `permeability` of `4π`.

$(textunits(Electronic,:Electronic))
""" Electronic

@doc """
    Natural::UnitSystem{1,1,1,1,1}

Natural `UnitSystem` with all primary constants having unit value.

$(textunits(Natural,:Natural))
""" Natural

@doc """
    NaturalGauss::UnitSystem{1,1,1,4π,1}

Natural (Gauss) `UnitSystem` with the Gaussian `permeability` value of `4π`.

$(textunits(NaturalGauss,:NaturalGauss))
""" NaturalGauss

@doc """
    QCD::UnitSystem{1,1,1,1,1/μₚₑ}

Qunatum chromodynamics `UnitSystem` with `electronmass` of `1/μₚₑ` or `1/$μₚₑ`.

$(textunits(QCD,:QCD))
""" QCD

@doc """
    QCDGauss::UnitSystem{1,1,1,4π,1/μₚₑ}

Qunatum chromodynamics (Gauss) `UnitSystem` with `electronmass` of `1/μₚₑ`.

$(textunits(QCDGauss,:QCDGauss))
""" QCDGauss

@doc """
    QCDoriginal::UnitSystem{1,1,1,4π/αinv,1/μₚₑ}

Qunatum chromodynamics (original) `UnitSystem` with `permeability` of `4π/αinv`.

$(textunits(QCDoriginal,:QCDoriginal))
""" QCDoriginal

# common conversions

@doc """
    kilograms(m::Real) = $(slug)m

Converts mass `m` from slugs to kilogram (kg).
""" kilograms, slug
@pure kilograms(m::Real,U::UnitSystem=English) = mass(m,Metric,U)

"""
    slugs(m::Real) = $(1/slug)m

Converts mass `m` from kilograms to slugs (slug).
"""
@pure slugs(m::Real,U::UnitSystem=Metric) = mass(m,English,U)

"""
    feet(d) = $(1/ft)d

Converts distance `d` from meters to feet (ft).
"""
@pure feet(d) = (lightspeed(English)/lightspeed(Metric))d

@doc """
    meters(d) = $(ft)d

Converts distance `d` from feet to meters (m).
""" meters, ft
@pure meters(d) = (lightspeed(Metric)/lightspeed(English))d

@doc """
    rankine*T = (9/5)*T

Converts temperature `T` from Kelvin to degrees Rankine (°R).
""" rankine

@doc """
    kelvin*T = (5/9)*T

Converts temperature `T` from degrees Rankine to Kelvin (K).
""" kelvin

"""
    moles(N::Real,U::UnitSystem=Metric) = N/avogadro(U)

Converts the number of molecules `N` to number of moles (mol).
"""
@pure moles(N::Real,U::UnitSystem) = N/avogadro(U)

"""
    molecules(n::Real,U::UnitSystem=Metric) = n*avogadro(U)

Converts the number of moles `n` to number of molecules (dimensionless).
"""
@pure molecules(n::Real,U::UnitSystem=Metric) = n*avogadro(U)

#const English = UnitSystem{49720.072683*mₑ/μₑᵤ/slug,ħ/slug/ft^2,𝘤/ft,4π,mₑ/slug}()
##const English = UnitSystem{rankine*1000/slug/ft*Rᵤ*mₑ/μₑᵤ,ħ/slug/ft^2,𝘤/ft,4π,mₑ/slug}()
#const English = UnitSystem{rankine*kB/ft/slug,ħ/slug/ft^2,𝘤/ft,4π,mₑ/slug}()
const EnglishNew = UnitSystem{universal(English)*boltzmann(English)*ft/rankine/Rᵤ,ħ/slug/ft^2,𝘤/ft,4π,mₑ/slug}()
@pure molarmass(U::UnitSystem{boltzmann(EnglishNew)}) = rankine*Rᵤ*mₑ/μₑᵤ/boltzmann(English)/slug/ft
