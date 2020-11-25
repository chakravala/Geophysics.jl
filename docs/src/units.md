# Dimensional unit systems

```@contents
Pages = ["units.md","references.md"]
```

A `UnitSystem` is a consistent set of dimensional values selected to accomodate a particular use case or standardization.
In total, five fundamental constants `kB,ħ,𝘤,μ,mₑ` are used to specify a specific unit system.
These are the constants of `boltzmann`, `planckreduced`, `lightspeed`, `permeability`, and `electronmass`.
Different choices of natural units or physical measurements result in a variety of unit systems optimized for many purposes.

```@docs
UnitSystem
```

## Metric SI Units

In the Systeme International d'Unites (the SI units) the `UnitSystem` constants are derived from the most accurate possible physical measurements and a few exactly defined constants.
Exact values are the `avogadro` number, `boltzmann` constant, `planck` constant, `lightspeed` definition, and elementary `charge` definition.

```Julia
julia> NA # avogadro
6.02214076e23

julia> kB # boltzmann
1.380649e-23

julia> 𝘩 # planck
6.62607015e-34

julia> 𝘤 # lightspeed
2.99792458e8

julia> 𝘦 # charge
1.602176634e-19
```

Physical measured values with uncertainty are the electron to proton mass ratio `μₑₐ`, proton to atomic mass ratio `μₚₐ`, fine structure constant `αinv`, the Rydberg `R∞` constant, and the Planck mass `mP`.

```Julia
julia> μₑₐ
0.0005485799090649074

julia> μₚₐ
1.007276466621

julia> αinv
137.035999084

julia> R∞ # rydberg
1.09737315681601e7

julia> mP # planckmass
2.176434e-8
```

From these numbers along with the `4π*1e-7` value of the Gaussian unit `μ₀`, the constants `planckreduced`, `permeability`, `electronmass`, `molarmass`, and proton to electon mass ratio are computed.

```Julia
julia> ħ # 𝘩/2π
1.0545718176461565e-34

julia> μ₀ # 2𝘩/𝘤/αinv/𝘦^2
1.256637062121048e-6

julia> mₑ # αinv^2*R∞*2𝘩/𝘤
9.109383701558256e-31

julia> Mᵤ # mₑ*NA/μₑₐ
0.000999999999656256

julia> μₚₑ # μₚₐ/μₑₐ
1836.152673432705
```

These result in variants based on the original `molarmass` constant and Gaussian `permeability` along with the 2019 redefined exact values.
```Julia
Metric::UnitSystem{Rᵤ*mₑ/μₑᵤ/0.001,ħ,𝘤,4π*1e-7,mₑ}
SI2019::UnitSystem{kB,ħ,𝘤,μ₀,mₑ}
```
```@docs
Geophysics.Metric
Geophysics.SI2019
```

Additional reference values include the ground state hyperfine structure transition frequency of caesium-133 `ΔνCs` and luminous efficacy `Kcd` of monochromatic radiation of 540 THz.

```Julia
julia> ΔνCs
9.19263177e9

julia> Kcd
683.0
```

## Other historic systems

Alternatives to the SI unit system are the centimetre-gram-second variants.
```Julia
CGS     ::UnitSystem{1e7*kB,1e7*ħ,100𝘤,4π,1000mₑ}
CGS2019 ::UnitSystem{1e7*kB,1e7*ħ,100𝘤,1e7*μ₀,1000mₑ}
```
```@docs
Geophysics.CGS
Geophysics.CGS2019
```

Historically, the `josephson` and `klitzing` constants have been used to define `Conventional` and `CODATA` derived `UnitSystem` variants.

```Julia
CODATA::UnitSystem{Rᵤ*mₑ/μₑᵤ/0.001,2/KJ2014/RK2014^2/π,𝘤,2KJ2014/𝘤/αinv,mₑ}()
Conventional::UnitSystem{1000Rᵤ*mₑ/μₑᵤ,2/KJ1990/RK1990^2/π,𝘤,2KJ1990/𝘤/αinv,mₑ}()
```

```Julia
julia> josephson(Conventional) # KJ1990
4.835979e14

julia> klitzing(Conventional) # RK1990
25812.807

julia> josephson(CODATA) # KJ2014
4.835978525e14

julia> klitzing(CODATA) # RK2014
25812.8074555
```

```@docs
Geophysics.Conventional
Geophysics.CODATA
```

In Britain and the United States an `English` system of engineering units was commonly used.

```@docs
Geophysics.English
```

## Natural units

With the introduction of the `planckmass` a set of natural atomic unit systems can be derived in terms of the gravitational coupling constant.

```Julia
julia> αG # (mₑ/mP)^2
1.751809945750515e-45
```

Some of the notable variants include

```Julia
Planck       ::UnitSystem{1,1,1,1,√(4π*αG)}
PlanckGauss  ::UnitSystem{1,1,1,4π,√αG}
Stoney       ::UnitSystem{1,αinv,1,4π,√(αG*αinv)}
Hartree      ::UnitSystem{1,1,αinv,4π/αinv^2,1}
Rydberg      ::UnitSystem{1,1,2αinv,π/αinv^2,1/2}
Schrodinger  ::UnitSystem{1,1,αinv,4π/αinv^2,√αinv*mₑ/mP}
Electronic   ::UnitSystem{1,αinv,1,4π,1}
Natural      ::UnitSystem{1,1,1,1,1}
NaturalGauss ::UnitSystem{1,1,1,4π,1}
QCD          ::UnitSystem{1,1,1,1,1/μₚₑ}
QCDGauss     ::UnitSystem{1,1,1,4π,1/μₚₑ}
QCDoriginal  ::UnitSystem{1,1,1,4π/αinv,1/μₚₑ}
```

```@docs
Geophysics.Planck
Geophysics.PlanckGauss
Geophysics.Stoney
Geophysics.Hartree
Geophysics.Rydberg
Geophysics.Schrodinger
Geophysics.Electronic
Geophysics.Natural
Geophysics.NaturalGauss
Geophysics.QCD
Geophysics.QCDGauss
Geophysics.QCDoriginal
```

## Fundamental constants of physics

The following are fundamental constants of physics:

```@docs
molarmass
avogadro
boltzmann
universal
planck
planckreduced
lightspeed
planckmass
atomicmass
protonmass
electronmass
newton
einstein
permeability
permittivity
coulomb
stefan
radiationintensity
impedance
charge
magneton
conductance
faraday
magneticflux
josephson
klitzing
hartree
rydberg
bohr
bohrreduced
electronradius
```

## Common conversion factors

Common conversion factors include `moles`, `molecules`, `kilograms`, `slugs`, `meters`, `feet`, `kelvin`, and `rankine`.

```@docs
moles
molecules
kilograms
slugs
meters
feet
kelvin
rankine
```

## Index

```@index
Pages = ["units.md"]
```

