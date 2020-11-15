#   This file is part of Geophysics.jl. It is licensed under the AGPL license
#   Geophysics Copyright (C) 2020 Michael Reed

export meters,feet, rankine, kelvin, moles, molecules, UnitSystem
export mass, slugs, kilograms

export CGS, CGS2019, Metric, SI2019, CODATA, Conventional, English
export Planck, PlanckGauss, Stoney, Hartree, Rydberg, Schrodinger, Electronic, Natural, NaturalGauss, QCD, QCDGauss, QCDoriginal

"""
    UnitSystem{kB,ħ,𝘤,μ₀,mₑ}

Standardized for engineering based on fundamental constants: `kB` Boltzmann's constant, `ħ` reduced Planck's constant, `𝘤` speed of light, `μ₀` vacuum permeability, and `mₑ` electron rest mass.
Primarily the `Metric` SI unit system is used in addition to the historic `English` engineering unit system.
These constants induce derived values for `avogadro`, `boltzmann`, `universal`, `planck`, `planckreduced`, `lightspeed`, `planckmass`, `atomicmass`, `protonmass`, `electronmass`, `newton`, `einstein`, `permeability`, `permittivity`, `coulomb`, and
additional constants `stefan`, `radiationintensity`, `impedance`, `charge`, `magneton`, `conductance`, `faraday`, `magneticflux`, `josephson`, `klitzing`, `hardtree`, `rydberg`, `bohr`, and `bohrreduced`.

Additional reference `UnitSystem` variants `CGS`, `CGS2019`, `SI2019`, `CODATA`, `Conventional`; along with several natural atomic units based on the fine structure constant `1/αinv` and the gravitational coupling constant `αG` (`Planck`, `PlanckGauss`, `Stoney`, `Hartree`, `Rydberg`, `Schrodinger`, `Electronic`, `Natural`, `NaturalGauss`, `QCD`, `QCDGauss`, and `QCDoriginal`).
""" #`Rᵤ,mᵤ,σ,ħ,μ₀,ε₀,kₑ,𝘦,𝔉,RK,Z₀,G₀`
struct UnitSystem{kB,ħ,𝘤,μ,mₑ} end
@pure boltzmann(::UnitSystem{k}) where k = k
@pure planckreduced(::UnitSystem{k,h}) where {k,h} = h
@pure lightspeed(::UnitSystem{k,h,c}) where {k,h,c} = c
@pure permeability(::UnitSystem{k,h,c,μ}) where {k,h,c,μ} = μ
@pure electronmass(::UnitSystem{k,h,c,μ,m}) where {k,h,c,μ,m} = m

@pure mass(U::UnitSystem,S::UnitSystem=Metric) = electronmass(U)/electronmass(S)
@pure mass(m::Real,U::UnitSystem=English,S::UnitSystem=Metric) = m*mass(U,S)
@pure planckmass(U::UnitSystem) = mass(mP,U)
@pure newton(U::UnitSystem) = lightspeed(U)*planckreduced(U)/planckmass(U)^2

Base.display(U::UnitSystem) = println("UnitSystem{kB=$(boltzmann(U)),ħ=$(planckreduced(U)),𝘤=$(lightspeed(U)),μ₀=$(permeability(U)),mᵤ=$(electronmass(U))}")

# fundamental constants

const NA,kB,𝘩,𝘤,𝘦 = 6.02214076e23,1.380649e-23,6.62607015e-34,299792458.,1.602176634e-19
const μₑₐ,μₚₐ,αinv,R∞ = 1/1822.888486209,1.007276466621,137.035999084,10973731.5681601
const Mᵤ = αinv^2*R∞*NA*2𝘩/𝘤/μₑₐ # Rydberg molar mass
const ħ,μ₀,mₑ,μₚₑ = 𝘩/2π,4π*1e-7,Mᵤ*μₑₐ/NA,μₚₐ/μₑₐ
const ΔνCs,Kcd,mP = 9192631770.0,683.0,2.176434e-8
const KJ1990,KJ2014 = 25812.807,25812.8074555
const δμ₀ = 2𝘩/𝘤/αinv/𝘦^2-μ₀ # ≈ 4π*5.5e-17, exact charge
const ħE = 0.7375598625957106ħ # × ft⋅lb⋅J⁻¹, mass(English)*feet(1)^2

# engineering units

const CGS = UnitSystem{1e7*kB,1e7*ħ,100𝘤,1e7*μ₀,1000mₑ}()
const CGS2019 = UnitSystem{1e7*kB,1e7*ħ,100𝘤,1e7(μ₀+δμ₀),1000mₑ}()
const Metric = UnitSystem{kB,ħ,𝘤,μ₀,mₑ}()
const SI2019 = UnitSystem{kB,ħ,𝘤,μ₀+δμ₀,mₑ}()
const CODATA = UnitSystem{kB,2/KJ2014/4.835978525e14^2/π,𝘤,2KJ2014/𝘤/αinv,mₑ}()
const Conventional = UnitSystem{kB,2/KJ1990/4.835979e14^2/π,𝘤,2KJ1990/𝘤/αinv,mₑ}()
const English = UnitSystem{5.657302466e-24,ħE,983571056.,μ₀,mₑ/0.45359237/32.17404856}()

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

@doc """
    planckmass(U::UnitSystem) = mass($(planckmass(Metric)),U)

Planck mass factor `mP` from the gravitational coupling constant `αG` (kg or slugs).
```Julia
juila> planckmass(Metric) # m³⋅kg⁻¹⋅s⁻²
$(planckmass(Metric))

julia> planckmass(English) # ft³⋅slug⁻¹⋅s⁻²
$(planckmass(English))
```
""" planckmass, mP

@doc """
    newton(x) = lightspeed(x)*planckreduced(x)/planckmass(x)^2

Universal gravitational constant `GG` of Newton's law (m³⋅kg⁻¹⋅s⁻² or ft³⋅slug⁻¹⋅s⁻²).
```Julia
juila> newton(Metric) # m³⋅kg⁻¹⋅s⁻²
$(newton(Metric))

julia> newton(English) # ft³⋅slug⁻¹⋅s⁻²
$(newton(English))
```
""" newton, GG

@doc """
    boltzmann(x) = universal(x)/avogadro(x)

Boltzmann constant `kB` is the entropy amount of a unit number microstate permutation.
```Julia
pressure*molecularmass == density*boltzmann*temperature
```
It satisfies the ideal gas law.

```Julia
julia> boltzmann(Metric) # J⋅K⁻¹
$(boltzmann(Metric))

julia> boltzmann(English) # ft⋅lb⋅°R⁻¹
$(boltzmann(English))
```
""" boltzmann, kB

@doc """
    planckreduced(x) = planck(x)/2π

Reduced Planck constant `ħ` is a Planck per radian (J⋅s⋅rad⁻¹ or ft⋅lb⋅s⋅rad⁻¹).

```Julia
julia> planckreduced(Metric) # J⋅s⋅rad⁻¹
$(planckreduced(Metric))

julia> planckreduced(English) # ft⋅lb⋅s⋅rad⁻¹
$(planckreduced(English))
```
""" planckreduced, ħ

@doc """
    lightspeed(x) = 1/sqrt(μ₀*ε₀)

Speed of light in a vacuum `𝘤` for massless particles (m⋅s⁻¹ or ft⋅s⁻¹).

```Julia
julia> universal(Metric) # m⋅s⁻¹
$(lightspeed(Metric))

julia> universal(English) # ft⋅s⁻¹
$(lightspeed(English))
```
""" lightspeed, 𝘤, cc

@doc """
    permeability(x) = 4π*1e-7

Magnetic permeability in a classical vacuum defined as `μ₀` in SI units (H⋅m⁻¹).

```Julia
julia> permeability(Metric) # H⋅m⁻¹
$(permeability(Metric))

julia> permeability(English) # slug⋅ft²⋅?⁻²
$(permeability(Metric))
```
""" permeability, μ₀, m0

@doc """
    electronmass(U::UnitSystem) = protonmass(U)/$μₚₑ

Electron rest mass unit `mₑ` of subatomic particle with `-𝘦` elementary charge  (kg or slugs).
```Julia
julia> electronmass(Metric) # kg
$(electronmass(Metric))

julia> electronmass(English) # slugs
$(electronmass(English))
```
""" electronmass, mₑ, me

@doc """
    Metric

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
    English

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
    CGS

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
    CGS2019

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
    CODATA

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
    SI2019

Systeme International d'Unites (the SI units) with `μ₀+$δμ₀` for a tuned `charge`.

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
""" SI2019

@doc """
    Conventional

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
    Planck

Planck `UnitSystem` with the `electronmass` value `√(4π*αG)` using gravitational coupling.

$(textunits(Planck,:Planck))
""" Planck

@doc """
    PlanckGauss

Planck (Gauss) `UnitSystem` with `permeability` of `4π` and `electronmass` coupling `√αG`.

$(textunits(PlanckGauss,:PlanckGauss))
""" PlanckGauss

@doc """
    Stoney

Stoney `UnitSystem` with `permeability` of `4π` and `electronmass` coupling `√(αG*αinv)`.

$(textunits(Stoney,:Stoney))
""" Stoney

@doc """
    Hartree

Hartree atomic `UnitSystem` with `lightspeed` of `αinv` and `permeability` of `4π/αinv^2`.

$(textunits(Hartree,:Hartree))
""" Hartree

@doc """
    Rydberg

Rydberg `UnitSystem` with `lightspeed` of `2αinv` and `permeability` of `π/αinv^2`.

$(textunits(Rydberg,:Rydberg))
""" Rydberg

@doc """
    Schrodinger

Schrodinger `UnitSystem` with `permeability` of `4π/αinv^2` and `electronmass` of `√αinv*mₑ/mP`.

$(textunits(Schrodinger,:Schrodinger))
""" Schrodinger

const Electronic = UnitSystem{1,αinv,1,4π,1}()
@doc """
    Electronic

Electronic `UnitSystem` with `planckreduced` of `αinv` and `permeability` of `4π`.

$(textunits(Electronic,:Electronic))
""" Electronic

@doc """
    Natural

Natural `UnitSystem` with all primary constants having unit value.

$(textunits(Natural,:Natural))
""" Natural

@doc """
    NaturalGauss

Natural (Gauss) `UnitSystem` with the Gaussian `permeability` value of `4π`.

$(textunits(NaturalGauss,:NaturalGauss))
""" NaturalGauss

@doc """
    QCD

Qunatum chromodynamics `UnitSystem` with `electronmass` of `1/μₚₑ` or `1/$μₚₑ`.

$(textunits(QCD,:QCD))
""" QCD

@doc """
    QCDGauss

Qunatum chromodynamics (Gauss) `UnitSystem` with `electronmass` of `1/μₚₑ`.

$(textunits(QCDGauss,:QCDGauss))
""" QCDGauss

@doc """
    QCDoriginal

Qunatum chromodynamics (original) `UnitSystem` with `permeability` of `4π/αinv`.

$(textunits(QCDoriginal,:QCDoriginal))
""" QCDoriginal

# constants

@pure avogadro(U::UnitSystem) = Mᵤ*μₑₐ/electronmass(U)
@doc """
    avogadro(x) = universal(x)/boltzmann(x) # Mᵤ/atomicmass(x), Mᵤ ≈ 0.001-3.5e-13

Avogadro `NA` is `molarmass(x)/molecularmass(x)` number of atoms in a 12 g sample of C₁₂.
```Julia
julia> avogadro(Metric) # mol⁻¹
$(avogadro(Metric))

julia> avogadro(English) # slug-mol⁻¹
$(avogadro(English))
```
""" avogadro, NA

@pure universal(U::UnitSystem) = boltzmann(U)*avogadro(U)
@doc """
    universal(x) = boltzmann(x)*avogadro(x)

Universal gas constant `Rᵤ` is factored into specific `gasconstant(x)*molarmass(x)` values.
```Julia
pressure*molarmass == density*universal*temperature
```
It satisfies the ideal gas law.

```Julia
julia> universal(Metric) # J⋅K⁻¹⋅mol⁻¹
$(universal(Metric))

julia> universal(English) # ft⋅lb⋅°R⁻¹⋅slug-mol⁻¹
$(universal(English))
```
""" universal, Rᵤ, Ru

@pure permittivity(U::UnitSystem) = inv(permeability(U)*lightspeed(U)^2)
@doc """
    permittivity(x) = 1/μ₀/lightspeed(x)^2

Dielectric permittivity constant `ε₀` of a classical vacuum (C²⋅N⁻¹⋅m⁻²).

```Julia
julia> permittivity(Metric) # C²⋅N⁻¹⋅m⁻²
$(permittivity(Metric))
```
""" permittivity, ε₀, e0

@pure coulomb(U::UnitSystem) = inv(4π*permittivity(U))
@doc """
    coulomb(x) = 1/4π/ϵ₀

Electrostatic proportionality constant `kₑ` for the Coulomb's law force (N⋅m²⋅C⁻²).

```Julia
julia> coulomb(Metric) # N⋅m²⋅C⁻²
$(coulomb(Metric))
```
""" coulomb, kₑ, ke

@pure planck(U::UnitSystem) = 2π*planckreduced(U)
@doc """
    planck(x) = 2π*planckreduced(x)

Planck constant `𝘩` is energy per electromagnetic frequency (J⋅s or ft⋅lb⋅s).

```Julia
julia> planck(Metric) # J⋅s
$(planck(Metric))

julia> planck(English) # ft⋅lb⋅s
$(planck(English))
```
""" planck, 𝘩, hh

@pure atomicmass(U::UnitSystem) = electronmass(U)/μₑₐ
@doc """
    atomicmass(U::UnitSystem) = Mᵤ/avogadro(U) # $Mᵤ ≈ 0.001-3.5e-13

Atomic mass unit `mᵤ` of 1/12 of the C₁₂ carbon-12 atom's mass  (kg or slugs).
```Julia
julia> atomicmass(Metric) # kg
$(atomicmass(Metric))

julia> atomicmass(English) # slugs
$(atomicmass(English))
```
""" atomicmass, mᵤ, mu

@pure protonmass(U::UnitSystem) =  μₚₐ*atomicmass(U)
@doc """
    protonmass(U::UnitSystem) = $(μₚₐ)atomicmass(U)

Proton mass unit `mₚ` of subatomic particle with `+𝘦` elementary charge  (kg or slugs).
```Julia
julia> protonmass(Metric) # kg
$(protonmass(Metric))

julia> protonmass(English) # slugs
$(protonmass(English))
```
""" protonmass, mₚ, mp

@pure einstein(U::UnitSystem) = 8π*newton(U)/lightspeed(U)^4
@doc """
    einstein(U::UnitSystem) = 8π*newton(U)/lightspeed(U)^4

Einstein's gravitational constant from the Einstein field equations (? or ?).
```Julia
julia> einstein(Metric) # ?
$(einstein(Metric))

julia> einstein(English) # ?
$(einstein(English))
```
""" einstein, κ

@pure stefan(U::UnitSystem) = 2π^5*boltzmann(U)^4/(15planck(U)^3*lightspeed(U)^2)
@doc """
    stefan(x) = 2π^5*boltzmann(U)^4/(15planck(U)^3*lightspeed(U)^2)

Stefan-Boltzmann proportionality `σ` of black body radiation (W⋅m⁻²⋅K⁻⁴ or ?⋅ft⁻²⋅°R⁻⁴).

```Julia
julia> stefan(Metric) # W⋅m⁻²⋅K⁻⁴
$(stefan(Metric))

julia> stefan(English) # lb⋅s⁻¹⋅ft⁻³⋅°R⁻⁴
$(stefan(English))
```
""" stefan, σ, SB

"""
    radiationdensity(x) = 4stefan(U)/lightspeed(U)

Raditation density constant of black body radiation (J⋅m⁻³⋅K⁻⁴ or lb⋅ft⁻²⋅°R⁻⁴).

```Julia
julia> radiationdensity(Metric) # J⋅m⁻³⋅K⁻⁴
$(radiationdensity(Metric))

julia> radiationdensity(English) # lb⋅ft⁻²⋅°R⁻⁴
$(radiationdensity(English))
```
"""
@pure radiationdensity(U::UnitSystem) = 4stefan(U)/lightspeed(U)

@pure impedance(U::UnitSystem) = permeability(U)*lightspeed(U)
@doc """
    impedance(U::UnitSystem) = permeability(U)*lightspeed(U)

Vacuum impedance of free space `Z₀` is magnitude ratio of electric to magnetic field (Ω).
```Julia
julia> impedance(Metric) # Ω
$(impedance(Metric))
```
""" impedance, Z₀, Z0

@pure charge(U::UnitSystem) = sqrt(2planck(U)/impedance(U)/αinv) # fine structure
@doc """
    charge(U::UnitSystem) = sqrt(2𝘩/$(αinv)impedance(U))

Quantized elementary charge `𝘦` of a proton or electron  (C).
```Julia
julia> charge(Metric) # C
$(charge(Metric))
```
""" charge, 𝘦, ee

@pure faraday(U::UnitSystem) = charge(U)*avogadro(U)
@doc """
    faraday(U::UnitSystem) = charge(U)*avogadro(U)

Electric charge per mole of electrons `𝔉` based on elementary charge (C⋅mol⁻¹).
```Julia
julia> faraday(Metric) # C⋅mol⁻¹
$(faraday(Metric))
```
""" faraday, 𝔉, FF

@pure josephson(U::UnitSystem) = 2charge(U)/planck(U)
@doc """
    josephson(U::UnitSystem) = 2charge(U)/planck(U)

Josephson constant `KJ` relating potential difference to irradiation frequency (Hz⋅V⁻¹).
```Julia
julia> josephson(Metric) # Hz⋅V⁻¹
$(josephson(Metric))
```
""" josephson, KJ

@pure magneticflux(U::UnitSystem) = inv(josephson(U))
@doc """
    magneticflux(U::UnitSystem) = planck(U)/2charge(U)

Magnetic flux quantum `Φ₀` is `1/josephson(U)` (Wb).
```Julia
julia> magneticflux(Metric) # Wb
$(magneticflux(Metric))
```
""" magneticflux, Φ₀

@pure klitzing(U::UnitSystem) = planck(U)/charge(U)^2
@doc """
    klitzing(U::UnitSystem) = 2/conductance(U)

Quantized Hall resistance `RK` (Ω).
```Julia
julia> klitzing(Metric) # Ω
$(klitzing(Metric))
```
""" klitzing, RK

@pure hardtree(U::UnitSystem) = electronmass(U)*(lightspeed(U)/αinv)^2
@doc """
    hardtree(U::UnitSystem) = electronmass(U)*(lightspeed(U)/$αinv)^2

Hardtree electric potential energy `Eₕ` of the hydrogen atom at ground state (J).
```Julia
julia> hardtree(Metric) # J
$(hardtree(Metric))
```
""" hardtree, Eₕ, Eh

@pure rydberg(U::UnitSystem) = hardtree(U)/2planck(U)/lightspeed(U)
@doc """
    rydberg(U::UnitSystem) = hardtree(U)/2planck(U)/lightspeed(U)

Rydberg constant `R∞` is lowest energy photon capable of ionizing atom at ground state (m⁻¹).
```Julia
julia> rydberg(Metric) # m⁻¹
$(rydberg(Metric))
```
""" rydberg, R∞

@pure bohr(U::UnitSystem) = αinv*planckreduced(U)/electronmass(U)/lightspeed(U)
@doc """
    bohr(U) = $αinv*planckreduced(U)/electronmass(U)/lightspeed(U)

Bohr radius of the hydrogen atom in its ground state `a₀` (m).
```Julia
julia> bohr(Metric) # m
$(bohr(Metric))
```
""" bohr, a₀, a0

"""
    bohrreduced(U::UnitSystem) = electronmass(U)/bohr(U)/$μₚₑ

Reduced Bohr radius including the effect of reduced mass in hydrogen atom (m).
```Julia
julia> bohrreduced(Metric) # m
$(bohrreduced(Metric))
```
"""
@pure bohrreduced(U::UnitSystem) = electronmass(U)*bohr(U)/μₚₑ

@pure electronradius(U::UnitSystem) = planckreduced(U)/electronmass(U)/lightspeed(U)/αinv
@doc """
    electronradius(U) = planckreduced(U)/electronmass(U)/lightspeed(U)/$αinv

Classical electron radius or Lorentz radius or Thomson scattering length (m).
```Julia
julia> electronradius(Metric) # m
$(electronradius(Metric))
```
""" electronradius, rₑ, re

@pure conductance(U::UnitSystem) = 2charge(U)^2/planck(U)
@doc """
    conductance(U::UnitSystem) = 2charge(U)^2/𝘩

Conductance quantum `G₀` is a quantized unit of electrical conductance (S).
```Julia
julia> conductance(Metric) # S
$(conductance(Metric))
```
""" conductance, G₀, G0

@pure magneton(U::UnitSystem) = charge(U)*planckreduced(U)/2electronmass(U)
"""
    magneton(U::UnitSystem) = charge(U)*planckreduced(U)/2electronmass(U)

Bohr magneton `μB` natural unit for expressing magnetic moment of electron (J⋅T⁻¹).
```Julia
julia> magneton(Metric) # J⋅T⁻¹
$(magneton(Metric))
```
""" magneton, μB

const κ = einstein(Metric)
const GG = newton(Metric)
const Rᵤ = universal(Metric)
const σ = stefan(Metric)
const μB = magneton(Metric)
const ε₀ = permittivity(Metric)
const kₑ = coulomb(Metric)
const mₚ = protonmass(Metric)
const mᵤ = atomicmass(Metric)
const 𝔉 = faraday(Metric)
const Φ₀ = magneticflux(Metric)
const Z₀ = impedance(Metric)
const G₀ = conductance(Metric)
const Eₕ = hardtree(Metric)
const a₀ = bohr(Metric)
const rₑ = electronradius(Metric)
const RK = klitzing(Metric)
const KJ = josephson(Metric)
const Mu,Ru,SB,hh,cc,m0,e0,ke,me,mp,mu,ee,FF,Z0,G0,Eh,a0,re = Mᵤ,Rᵤ,σ,𝘩,𝘤,μ₀,ε₀,kₑ,mₑ,mₚ,mᵤ,𝘦,𝔉,Z₀,G₀,Eₕ,a₀,rₑ
export κ, GG, NA, kB, Rᵤ, σ, 𝘩, ħ, 𝘤, μ₀, ε₀, kₑ, mₑ, mₚ, mᵤ, 𝘦, 𝔉, Φ₀, Z₀, G₀, Eₕ, R∞, a₀, rₑ, KJ, RK, Ru, SB, hh, cc, m0, e0, ke, me, mp, mu, ee, FF, Z0, G0, Eh, a0, re
export αG, αinv, μₚₑ, μₑₐ, μₚₐ, mpe, mea, mpa, mP, δμ₀, Mᵤ, Mu
const mpe, mea, mpa, SI = μₚₑ, μₑₐ, μₚₐ, ΔνCs, Kcd, SI2019

export electronmass, protonmass, atomicmass, planckmass, stefan, radiationintensity, einstein, impedance, charge, faraday, josephson, klitzing, hardtree, rydberg, bohr, bohrreduced, electronradius, conductance, magneticflux, magneton

const Constants = (:newton,:avogadro,:boltzmann,:planck,:planckreduced,:lightspeed,:universal,:permeability,:permittivity,:coulomb)

const Properties = (:units,:molarmass,:molecularmass,:gasconstant,Constants...)

const Intrinsic = (:viscosity,:conductivity,:heatvolume,:heatpressure,:heatratio,:prandtl,:sonicspeed,:freedom,:energy,:enthalpy)

"""
    kilograms(m::Real) = $(kilograms(1))m

Converts mass `m` from slugs to kilogram (kg).
"""
@pure kilograms(m::Real,U::UnitSystem=English) = mass(m,Metric,U)

"""
    slugs(m::Real) = $(slugs(1))m

Converts mass `m` from kilograms to slugs (slug).
"""
@pure slugs(m::Real,U::UnitSystem=Metric) = mass(m,English,U)

"""
    feet(d) = $(feet(1))d

Converts distance `d` from meters to feet (ft).
"""
@pure feet(d) = (lightspeed(English)/lightspeed(Metric))d

"""
    meters(d) = $(meters(1))d

Converts distance `d` from feet to meters (m).
"""
@pure meters(d) = (lightspeed(Metric)/lightspeed(English))d

"""
    rankine(T) = (9/5)T

Converts temperature `T` from Kelvin to degrees Rankine (°R).
"""
@pure rankine(T) = (9/5)T

"""
    kelvin(T) = (5/9)T

Converts temperature `T` from degrees Rankine to Kelvin (K).
"""
@pure kelvin(T) = (5/9)T

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
