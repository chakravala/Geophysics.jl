# The UnitSystem

*Physical unit system constants (Metric, English, Natural, etc...)*

[![DOI](https://zenodo.org/badge/317419353.svg)](https://zenodo.org/badge/latestdoi/317419353)

```@contents
Pages = ["units.md","constants.md","convert.md"]
```

Specifications for dimensional units are in the [UnitSystems.jl](https://github.com/chakravala/UnitSystems.jl) and [UnitfulSystems.jl](https://github.com/chakravala/UnitfulSystems.jl) repositories.
The two packages are designed so that they can be interchanged if compatibility with [Unitful.jl](https://github.com/PainterQubits/Unitful.jl) is desired or not.
However, the `UnitfulSystems` package has fewer `UnitSystem` specifications available than the `UnitSystems` package due to limitations in combination with the `Unitful` package.
Specifically, `Metric`, `SI2019`, `CODATA`, `Conventional`, `MTS`, `EMU2019`, `English`, and `EnglishUS` can have `Unitful` values; while `Gauss`, `LorentzHeaviside`, `Thomson`, `EMU`, `ESU`, `ESU2019`, `IAU`, `FFF`, `Planck`, `PlanckGauss`, `Stoney`, `Hartree`, `Rydberg`, `Schrodinger`, `Electronic`, `Natural`, `NaturalGauss`, `QCD`, `QCDGauss`, and `QCDoriginal` currently only support plain numerical values.

```Julia
pkg> add UnitSystems # or UnitfulSystems

julia> using UnitSystems
```

A `UnitSystem` is a consistent set of dimensional values selected to accomodate a particular use case or standardization.
It is possible to convert derived physical quantities from any `UnitSystem` specification into any other using accurate values.
In total, five fundamental constants `kB,ħ,𝘤,μ₀,mₑ` are used to specify a specific unit system.
These are the constants of `boltzmann`, `planckreduced`, `lightspeed`, `permeability`, and `electronmass`.
Different choices of natural units or physical measurements result in a variety of unit systems optimized for many purposes.

```math
k_B, \qquad \hbar, \qquad c, \qquad \mu_0, \qquad m_e, \qquad (M_u), \qquad (\lambda), \qquad (\alpha_L)
```
Another important additional definition is the `molarmass` constant `Mᵤ`, which is automatically selected based on the choice of `boltzmann` constant (but can also be customized if necessary).
Historically, older electromagnetic unit systems also relied on a `rationalization` constant `λ` and a `lorentz` force proportionality constant `αL`.
In most unit systems these extra constants have a value of `1` unless otherwise specified.

```@docs
UnitSystem
```

Other similar packages include [PhysicalConstants.jl](https://github.com/JuliaPhysics/PhysicalConstants.jl), [MathPhysicalConstants.jl](https://github.com/LaGuer/MathPhysicalConstants.jl), [Unitful.jl](https://github.com/PainterQubits/Unitful.jl.git), [UnitfulSystems.jl](https://github.com/chakravala/UnitfulSystems.jl), [UnitfulUS.jl](https://github.com/PainterQubits/UnitfulUS.jl), [UnitfulAstro.jl](https://github.com/JuliaAstro/UnitfulAstro.jl), [UnitfulAtomic.jl](https://github.com/sostock/UnitfulAtomic.jl), [NaturallyUnitful.jl](https://github.com/MasonProtter/NaturallyUnitful.jl), and [UnitfulMoles.jl](https://github.com/rafaqz/UnitfulMoles.jl).

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

Physical measured values with uncertainty are electron to proton mass ratio `μₑᵤ`, proton to atomic mass ratio `μₚᵤ`, inverted fine structure constant `αinv`, the Rydberg `R∞` constant, and the Planck mass `mP`.

```math
\mu_{eu} = \frac{m_e}{m_u} \approx \frac{1}{1822.9},
\mu_{pu} = \frac{m_p}{m_u} \approx 1.00727647,
\alpha \approx \frac{1}{137.036},
R_\infty \approx 1.097373\mathrm{e}{7},
m_P \approx 2.176434\mathrm{e}{-8},
```

```Julia
julia> μₑᵤ # mₑ/mᵤ
0.0005485799090649074

julia> μₚᵤ # mₚ/mᵤ
1.007276466621

julia> αinv # 1/(fine structure)
137.035999084

julia> R∞ # rydberg
1.09737315681601e7

julia> mP # planckmass
2.176434e-8
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
julia> ħ # 𝘩/2π
1.0545718176461565e-34

julia> μ₀ # 2𝘩/𝘤/αinv/𝘦^2
1.256637062121048e-6

julia> mₑ # αinv^2*R∞*2𝘩/𝘤
9.109383701558256e-31

julia> Mᵤ # mₑ*NA/μₑᵤ
0.000999999999656256

julia> μₚₑ # μₚᵤ/μₑᵤ, mₚ/mₑ
1836.152673432705
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

```Julia
Metric::UnitSystem{Rᵤ*mₑ/μₑᵤ/0.001,ħ,𝘤,4π*1e-7,mₑ}
SI2019::UnitSystem{kB,ħ,𝘤,μ₀,mₑ}
```

```math
\frac{1}{1000} - M_u = 3.437439135417497\mathrm{e}{-13}, \qquad
\frac{4π}{10^7} - \mu_0 \approx -6.851306461996397\mathrm{e}{-16}
```

```@docs
Metric
SI2019
```

Additional reference values include the ground state `hyperfine` structure transition frequency of caesium-133 `ΔνCs` and `luminousefficacy` of monochromatic radiation `Kcd` of 540 THz.

```Julia
julia> ΔνCs # hyperfine
9.19263177e9

julia> Kcd # luminousefficacy
683.002
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
```Julia
EMU              ::UnitSystem{1e10*Rᵤ*mₑ/μₑᵤ,1e7*ħ,100𝘤,1,1000mₑ,4π}
ESU              ::UnitSystem{1e10*Rᵤ*mₑ/μₑᵤ,1e7*ħ,100𝘤,(100𝘤)^-2,1000mₑ,4π}
Gauss            ::UnitSystem{1e10*Rᵤ*mₑ/μₑᵤ,1e7*ħ,100𝘤,1,1000mₑ,4π,0.01/𝘤}
LorentzHeaviside ::UnitSystem{1e10*Rᵤ*mₑ/μₑᵤ,1e7*ħ,100𝘤,1,1000mₑ,1,0.01/𝘤}
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
```
Newer modern and rationalized variants of electromagnetic and electrostatic units are also made available.
```Julia
EMU2019::UnitSystem{1e7*kB,1e7*ħ,100𝘤,1e7*μ₀,1000mₑ}
ESU2019::UnitSystem{1e7*kB,1e7*ħ,100𝘤,1e3*μ₀/𝘤^2,1000mₑ}
```
```@docs
EMU2019
ESU2019
```

## Historical Unit Systems

Historically, the `josephson` and `klitzing` constants have been used to define `Conventional` and `CODATA` variants.

```math
\tilde k_B = \frac{8 \tilde R_\infty R_u}{c\alpha^2\mu_{eu} \tilde M_u \tilde K_J^2 \tilde R_K}, \quad
\tilde \hbar = \frac{2}{\pi\tilde K_J^2 \tilde R_K}, \quad
\tilde c = c, \quad
\tilde\mu_0 = \frac{2\tilde R_K\alpha}{c}, \quad
\tilde m_e = \frac{8\tilde R_\infty}{\tilde K_J^2\tilde R_Kc\alpha^2}, \quad
(\tilde M_u = \frac{1}{1000})
```

```Julia
CODATA       ::UnitSystem{Rᵤ*mₑ2014/μₑᵤ/0.001,2/RK2014/KJ2014^2/π,𝘤,2RK2014/𝘤/αinv,mₑ2014}
Conventional ::UnitSystem{Rᵤ*mₑ1990/μₑᵤ/0.001,2/RK1990/KJ1990^2/π,𝘤,2RK1990/𝘤/αinv,mₑ1990}
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
```

In Britain and the United States an `English` system of engineering units was commonly used.

```math
\tilde k_B = \frac{m_e R_u}{\mu_{eu} M_u\text{slug}\,\text{ft}^2}, \quad
\tilde \hbar = \frac{\hbar}{\text{slug}\cdot \text{ft}^2}, \quad
\tilde c = \frac{c}{\text{ft}}, \quad
\tilde\mu_0 = 4\pi, \quad
\tilde m_e = \frac{m_e}{\text{slug}}, \quad
(\tilde M_u = 1)
```
```Julia
English   ::UnitSystem{kB*rankine/slug/ft^2,ħ/slug/ft^2,𝘤/ft,4π,mₑ/slug}
EnglishUS ::UnitSystem{1000Rᵤ*mₑ/μₑᵤ*rankine/slug/ftUS^2,ħ/slug/ftUS^2,𝘤/ftUS,4π,mₑ/slug}
```
```@docs
English
EnglishUS
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
```

An impractical yet humorous unit system is the `FFF` specification.
```@docs
FFF
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
1.751809945750515e-45
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

