# The UnitSystem

*Physical unit system constants (Metric, English, Natural, etc...)*

[![DOI](https://zenodo.org/badge/317419353.svg)](https://zenodo.org/badge/latestdoi/317419353)

```@contents
Pages = ["units.md","constants.md","convert.md"]
```

> In fact there is nothing transcendental about dimensions; the ultimate principle is precisely expressible (in Newton's terminology) as one of *similitude*, exact or approximate, to be tested by the rule that mere change in the magnitudes of the ordered scheme of units of measurement that is employed must not affect sensibly the forms of the equations that are the adequate expression of the underlying relations of the problem. (J.L.)

Specifications for dimensional units are in the [UnitSystems.jl](https://github.com/chakravala/UnitSystems.jl) and [Similitude.jl](https://github.com/chakravala/Similitude.jl) and [MeasureSystems.jl](https://github.com/chakravala/MeasureSystems.jl) repositories.
The three packages are designed so that they can be interchanged with compatibility.
On its own `UnitSystems` is the fastest package, while `Similitude` (provides `Quantity` type) and `MeasureSystems` (introduces [Measurements.jl](https://github.com/JuliaPhysics/Measurements.jl) uncertainty) build additional features on top of `UnitSystems` base defintions.
Additionally, in the `UnitSystems` repository there is an equivalent [Wolfram language paclet](https://reference.wolfram.com/language/guide/Paclets) `Kernel` and also an unmaintained Rust `src` implementation.
Defaults are shared across the packages: `Metric`, `SI2019`, `CODATA`, `Conventional`, `International`, `InternationalMean`, `MetricEngineering`, `SI2019Engineering`, `GravitationalMetric`, `GravitationalSI2019`, `British`, `British2019`, `Survey`, `Survey2019`, `English`, `English2019`, `FPS`, `FPS2019`, `Gauss`, `LorentzHeaviside`, `Thomson`, `EMU`, `ESU`, `EMU2019`, `ESU2019`, `IAU`, `IAUE`, `IAUJ`, `Astronomical`, `Hubble`, `Cosmological`, `CosmologicalQuantum`, `Nautical`, `MPH`, `KKH`, `MTS`, `FFF`, `Planck`, `PlanckGauss`, `Stoney`, `Hartree`, `Rydberg`, `Schrodinger`, `Electronic`, `Natural`, `NaturalGauss`, `QCD`, `QCDGauss`, and `QCDoriginal`.

```Julia
julia> using UnitSystems # or Similitude or MeasureSystems
```

A `UnitSystem` is a consistent set of dimensional values selected to accomodate a particular use case standardization.
It is possible to convert derived physical quantities from any `UnitSystem` specification into any other using accurate values.
Eleven fundamental constants `kB`, `ħ`, `𝘤`, `μ₀`, `mₑ`, `Mᵤ`, `Kcd`, `θ`, `λ`, `αL`, `g₀` are used to govern a specific unit system consistent scaling.
These are the constants `boltzmann`, `planckreduced`, `lightspeed`, `vacuumpermeability`, `electronmass`, `molarmass`, `luminousefficacy`, `angle`, `rationalization`, `lorentz`, and `gravity`.
Different choices of natural units or physical measurements result in a variety of unit systems for many purposes.

```math
k_B, \qquad \hbar, \qquad c, \qquad \mu_0, \qquad m_e, \qquad M_u, \qquad K_{cd}, \qquad \theta, \qquad \lambda, \qquad \alpha_L, \qquad g_0
```
Historically, older electromagnetic unit systems also relied on a `rationalization` constant `λ` and a `lorentz` force proportionality constant `αL`.
In most unit systems these extra constants have a value of `1` unless otherwise specified.

```@docs
MeasureSystems.UnitSystem
```

Specification of `Universe` with the dimensionless `Coupling` constants `coupling`, `finestructure`, `electronunit`, `protonunit`, `protonelectron`, and `darkenergydensity`. Alterations to these values can be facilitated and quantified using parametric polymorphism.
Due to the `Coupling` interoperability, the `MeasureSystems` package is made possible to support calculations with `Measurements` having error standard deviations.

Similar packages: [UnitSystems.jl](https://github.com/chakravala/UnitSystems.jl), [Similitude.jl](https://github.com/chakravala/Similitude.jl), [MeasureSystems.jl](https://github.com/chakravala/MeasureSystems.jl), [PhysicalConstants.jl](https://github.com/JuliaPhysics/PhysicalConstants.jl), [MathPhysicalConstants.jl](https://github.com/LaGuer/MathPhysicalConstants.jl), [Unitful.jl](https://github.com/PainterQubits/Unitful.jl), [UnitfulUS.jl](https://github.com/PainterQubits/UnitfulUS.jl), [UnitfulAstro.jl](https://github.com/JuliaAstro/UnitfulAstro.jl), [UnitfulAtomic.jl](https://github.com/sostock/UnitfulAtomic.jl), [NaturallyUnitful.jl](https://github.com/MasonProtter/NaturallyUnitful.jl), and [UnitfulMoles.jl](https://github.com/rafaqz/UnitfulMoles.jl).

### Default UnitSystems

```@index
Pages = ["units.md"]
```

## Metric SI Unit Systems

In the Systeme International d'Unites (the SI units) the `UnitSystem` constants are derived from the most accurate possible physical measurements and a few exactly defined constants.
Exact values are the `avogadro` number, `boltzmann` constant, `planck` constant, `lightspeed` definition, and elementary `charge` definition.

```math
N_A = 6.02214076\mathrm{e}{23},
k_B = 1.380649\mathrm{e}{-23},
h = 6.62607015\mathrm{e}{-34},
c = 299792458,
e = 1.602176634\mathrm{e}{-19}
```

```Julia
julia> NA # avogadro
NA = 6.02214076e23

julia> kB # boltzmann
kB = 1.380649e-23

julia> 𝘩 # planck
𝘩 = 6.62607015e-34

julia> 𝘤 # lightspeed
𝘤 = 2.99792458e8

julia> 𝘦 # charge
𝘦 = 1.602176634e-19
```

Physical measured values with uncertainty are electron to proton mass ratio `μₑᵤ`, proton to atomic mass ratio `μₚᵤ`, inverted fine structure constant `αinv`, the Rydberg `R∞` constant, and the Planck mass `mP`.

```math
\mu_{eu} = \frac{m_e}{m_u} \approx \frac{1}{1822.9},
\mu_{pu} = \frac{m_p}{m_u} \approx 1.00727647,
\alpha \approx \frac{1}{137.036},
R_\infty \approx 1.097373\mathrm{e}{7},
m_P \approx 2.176434\mathrm{e}{-8},
```

```Julia
julia> μₑᵤ # electronunit
μₑᵤ = 0.000548579909065 ± 1.6e-14

julia> μₚᵤ # protonunit
μₚᵤ = 1.007276466621 ± 5.3e-11

julia> αinv # 1/finestructure
α⁻¹ = 137.035999084 ± 2.1e-8

julia> R∞ # rydbgerg
R∞ = 1.097373156816e7 ± 2.1e-5

julia> mP # planckmass
mP = 2.176434e-8 ± 2.4e-13
```

From these numbers along with the optional `4π` Gaussian `rationalization` value, the constants `planckreduced`, `permeability`, `electronmass`, `molarmass`, and proton to electon mass ratio are computed.

```math
\hbar = \frac{h}{2\pi}, \qquad
\mu_0 = \frac{2h\alpha}{ce^2}, \qquad
m_e = \frac{2hR_\infty}{c\alpha^2}, \qquad
M_u = \frac{m_e}{\mu_{eu}}N_A = \frac{2h R_\infty N_A}{c\alpha^2\mu_{eu}}, \qquad
\mu_{pe} = \frac{\mu_{pu}}{\mu_{eu}} = \frac{m_p}{m_e}
```

```Julia
julia> ħ # planckreduced
𝘩*τ⁻¹ = 1.0545718176461565e-34

julia> μ₀ # vacuumpermeability
𝘩*𝘤⁻¹𝘦⁻²α*2 = 1.25663706212e-6 ± 1.9e-16

julia> mₑ # electronmass
𝘩*𝘤⁻¹R∞*α⁻²2 = 9.1093837016e-31 ± 2.8e-40

julia> Mᵤ # molarmass
𝘩*𝘤⁻¹NA*R∞*α⁻²μₑᵤ⁻¹2 = 0.00099999999966 ± 3.1e-13

julia> μₚₑ # protonelectron
μₑᵤ⁻¹μₚᵤ = 1836.15267343 ± 1.1e-7
```

These result in variants based on the original `molarmass` constant and Gaussian `permeability` along with the 2019 redefined exact values. The main difference between the two is determined by $\delta\tilde M_u$ and $\delta\tilde\mu_0$.

```math
\tilde k_B = \frac{m_e R_u}{\mu_{eu} \tilde M_u} = \frac{k_B N_A}{m_u \tilde M_u}, \quad
\tilde \hbar = \hbar, \quad
\tilde c = c, \quad
\tilde\mu_0 = \frac{4\pi}{10^7} + \delta\tilde \mu_0, \quad
\tilde m_e = m_e, \quad
(\tilde M_u = \frac{1}{1000} + \delta \tilde M_u)
```

```@docs
MetricSystem
MeasureSystems.Metric
SI2019
MetricEngineering
SI2019Engineering
SI1976
```

Additional reference values include the ground state `hyperfine` structure transition frequency of caesium-133 `ΔνCs` and `luminousefficacy` of monochromatic radiation `Kcd` of 540 THz.

```Julia
julia> ΔνCs # hyperfine
ΔνCs = 9.19263177e9

julia> Kcd # luminousefficacy
Kcd = 683.01969009009
```

## Electromagnetic CGS Systems

Alternatives to the SI unit system are the centimetre-gram-second variants.

```math
\tilde k_B = 10^7\frac{\tilde m_e R_u}{\mu_{eu} \tilde M_u}, \quad
\tilde \hbar = 10^7\hbar, \quad
\tilde c = 100c, \quad
\tilde\mu_0 = 4\pi + \delta\tilde \mu_0, \quad
\tilde m_e = 1000m_e, \quad
(\tilde M_u, \, \tilde \lambda,\, \tilde \alpha_L)
```
There are multiple choices of elctromagnetic units for these variants based on electromagnetic units, electrostatic units, Gaussian non-rationalized units, and Lorentz-Heaviside rationalized units.
```@docs
GaussSystem
```
Note that `CGS` is an alias for the `Gauss` system.
```@docs
EMU
ESU
Gauss
LorentzHeaviside
```
When `Thomson` originally derived Maxwell's equations using electromagnetic notation, he arrived at a factor of `1/2` for the `lorentz` force constant, resulting in a slightly different sytem.
```@docs
Thomson
Kennelly
```
## Modified (Entropy) Unit Systems

```@docs
EntropySystem
GravitationalMetric
GravitationalSI2019
```
Newer modern and rationalized variants of electromagnetic and electrostatic units are also made available.
```@docs
EMU2019
ESU2019
ElectricSystem
International
InternationalMean
```

Historically, the `josephson` and `klitzing` constants have been used to define `Conventional` and `CODATA` variants.

```math
\tilde k_B = \frac{8 \tilde R_\infty R_u}{c\alpha^2\mu_{eu} \tilde M_u \tilde K_J^2 \tilde R_K}, \quad
\tilde \hbar = \frac{2}{\pi\tilde K_J^2 \tilde R_K}, \quad
\tilde c = c, \quad
\tilde\mu_0 = \frac{2\tilde R_K\alpha}{c}, \quad
\tilde m_e = \frac{8\tilde R_\infty}{\tilde K_J^2\tilde R_Kc\alpha^2}, \quad
(\tilde M_u = \frac{1}{1000})
```

```@docs
ConventionalSystem
```

```Julia
julia> josephson(Conventional) # KJ1990
KJ90 = 4.835979e14 [M⁻¹L⁻²TQ] Conventional

julia> klitzing(Conventional) # RK1990
RK90 = 25812.807 [ML²T⁻¹Q⁻²] Conventional

julia> josephson(CODATA) # KJ2014
KJ = 4.835978525e14 ± 3.0e6 [M⁻¹L⁻²TQ] CODATA

julia> klitzing(CODATA) # RK2014
RK = 25812.8074555 ± 5.9e-6 [ML²T⁻¹Q⁻²] CODATA
```

```@docs
Conventional
CODATA
```

In the Soviet Union, a metre-tonne-second system was also used briefly.
```math
\tilde k_B = 10^3\frac{\tilde m_e R_u}{\mu_{eu} \tilde M_u}, \quad
\tilde \hbar = 1000\hbar, \quad
\tilde c = c, \quad
\tilde\mu_0 = \frac{4\pi}{1000}, \quad
\tilde m_e = \frac{m_e}{1000}, \quad
(\tilde M_u = 10^{-6})
```
```@docs
MTS
KKH
MPH
Nautical
```

## Foot-Pound-Second-Rankine

In Britain and the United States an `English` system of engineering units was commonly used.

```math
\tilde k_B = \frac{m_e R_u}{\mu_{eu} M_u\text{slug}\,\text{ft}^2}, \quad
\tilde \hbar = \frac{\hbar}{\text{slug}\cdot \text{ft}^2}, \quad
\tilde c = \frac{c}{\text{ft}}, \quad
\tilde\mu_0 = 4\pi, \quad
\tilde m_e = \frac{m_e}{\text{slug}}, \quad
(\tilde M_u = 1)
```

```@docs
RankineSystem
MeasureSystems.British
British2019
Survey
Survey2019
MeasureSystems.English
MeasureSystems.English2019
FPS
FPS2019
```

An impractical yet humorous unit system is the `FFF` specification.
```@docs
FFF
```

## Astronomical Unit Systems

```@docs
Astronomical
```

The International Astronomical Union (IAU) units are based on the solar mass, distance from the sun to the earth, and the length of a terrestrial day.

```math
\tilde k_B = \frac{\tilde m_e R_u\text{day}^2}{\mu_{eu} \tilde M_u\text{au}^2}, \quad
\tilde \hbar = \frac{\hbar \text{day}}{m_\odot \text{au}^2}, \quad
\tilde c = c\frac{\text{day}}{\text{au}}, \quad
\tilde\mu_0 = \frac{4\pi}{10^7 m_\odot \text{au}^2}, \quad
\tilde m_e = \frac{m_e}{m_\odot}, \quad
(\tilde M_u = \frac{1}{1000m_\odot})
```

```@docs
IAU
IAUE
IAUJ
Hubble
Cosmological
CosmologicalQuantum
```

## Natural Unit Systems

With the introduction of the `planckmass` a set of natural atomic unit systems can be derived in terms of the gravitational coupling constant.

```math
\alpha_G = \left(\frac{m_e}{m_P}\right)^2, \qquad
\tilde k_B = 1, \qquad
(\tilde M_u = 1, \quad \tilde \lambda = 1, \quad \tilde \alpha_L = 1)
```

```Julia
julia> αG # (mₑ/mP)^2
𝘩²𝘤⁻²mP⁻²R∞²α⁻⁴2² = 1.75181e-45 ± 3.9e-50
```

Some of the notable variants include

```Julia
Planck       ::UnitSystem{1,1,1,1,√(4π*αG)}
PlanckGauss  ::UnitSystem{1,1,1,4π,√αG}
Stoney       ::UnitSystem{1,αinv,1,4π,√(αG*αinv)}
Hartree      ::UnitSystem{1,1,αinv,4π/αinv^2,1}
Rydberg      ::UnitSystem{1,1,2αinv,π/αinv^2,1/2}
Schrodinger  ::UnitSystem{1,1,αinv,4π/αinv^2,√(αG*αinv)}
Electronic   ::UnitSystem{1,αinv,1,4π,1}
Natural      ::UnitSystem{1,1,1,1,1}
NaturalGauss ::UnitSystem{1,1,1,4π,1}
QCD          ::UnitSystem{1,1,1,1,1/μₚₑ}
QCDGauss     ::UnitSystem{1,1,1,4π,1/μₚₑ}
QCDoriginal  ::UnitSystem{1,1,1,4π/αinv,1/μₚₑ}
```

```math
\tilde k_B = 1, \qquad
\tilde \hbar = 1, \qquad
\tilde c = 1, \qquad
\tilde \mu_0 = 1, \qquad
\tilde m_e = \sqrt{4\pi \alpha_G}
```
```@docs
Planck
```

```math
\tilde k_B = 1, \qquad
\tilde \hbar = 1, \qquad
\tilde c = 1, \qquad
\tilde \mu_0 = 4\pi, \qquad
\tilde m_e = \sqrt{4\pi \alpha_G}
```
```@docs
PlanckGauss
```

```math
\tilde k_B = 1, \qquad
\tilde \hbar = \frac1\alpha, \qquad
\tilde c = 1, \qquad
\tilde \mu_0 = 4\pi, \qquad
\tilde m_e = \sqrt{\frac{\alpha_G}{\alpha}}
```
```@docs
Stoney
```

```math
\tilde k_B = 1, \qquad
\tilde \hbar = 1, \qquad
\tilde c = \frac1\alpha, \qquad
\tilde \mu_0 = 4\pi\alpha^2, \qquad
\tilde m_e = 1
```
```@docs
Hartree
```

```math
\tilde k_B = 1, \qquad
\tilde \hbar = 1, \qquad
\tilde c = \frac2\alpha, \qquad
\tilde \mu_0 = \pi\alpha^2, \qquad
\tilde m_e = \frac{1}{2}
```
```@docs
Rydberg
```

```math
\tilde k_B = 1, \qquad
\tilde \hbar = 1, \qquad
\tilde c = \frac1\alpha, \qquad
\tilde \mu_0 = 4\pi\alpha^2, \qquad
\tilde m_e = \sqrt{\frac{\alpha_G}{\alpha}}
```
```@docs
Schrodinger
```

```math
\tilde k_B = 1, \qquad
\tilde \hbar = \frac1\alpha, \qquad
\tilde c = 1, \qquad
\tilde \mu_0 = 4\pi, \qquad
\tilde m_e = 1
```
```@docs
Electronic
```

```math
\tilde k_B = 1, \qquad
\tilde \hbar = 1, \qquad
\tilde c = 1, \qquad
\tilde \mu_0 = 1, \qquad
\tilde m_e = 1
```
```@docs
Natural
```

```math
\tilde k_B = 1, \qquad
\tilde \hbar = 1, \qquad
\tilde c = 1, \qquad
\tilde \mu_0 = 4\pi, \qquad
\tilde m_e = 1
```
```@docs
NaturalGauss
```

```math
\tilde k_B = 1, \qquad
\tilde \hbar = 1, \qquad
\tilde c = 1, \qquad
\tilde \mu_0 = 1, \qquad
\tilde m_e = \frac1{\mu_{pe}} = \frac{m_e}{m_p}
```
```@docs
QCD
```

```math
\tilde k_B = 1, \qquad
\tilde \hbar = 1, \qquad
\tilde c = 1, \qquad
\tilde \mu_0 = 4\pi, \qquad
\tilde m_e = \frac1{\mu_{pe}} = \frac{m_e}{m_p}
```
```@docs
QCDGauss
```

```math
\tilde k_B = 1, \qquad
\tilde \hbar = 1, \qquad
\tilde c = 1, \qquad
\tilde \mu_0 = 4\pi\alpha, \qquad
\tilde m_e = \frac1{\mu_{pe}} = \frac{m_e}{m_p}
```
```@docs
QCDoriginal
```

## UnitSystem Index

```@index
Pages = ["units.md","constants.md"]
```

