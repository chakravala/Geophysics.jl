# Dimensional unit systems

```@contents
Pages = ["units.md","references.md"]
```

Specifications for dimensional units are in the [UnitSystems.jl](https://github.com/chakravala/UnitSystems.jl) repository.
```Julia
pkg> add UnitSystems

julia> using UnitSystems
```

A `UnitSystem` is a consistent set of dimensional values selected to accomodate a particular use case or standardization.
In total, five fundamental constants `kB,ħ,𝘤,μ₀,mₑ` are used to specify a specific unit system.
These are the constants of `boltzmann`, `planckreduced`, `lightspeed`, `permeability`, and `electronmass`.
Different choices of natural units or physical measurements result in a variety of unit systems optimized for many purposes.

```math
k_B, \qquad \hbar, \qquad c, \qquad \mu_0, \qquad m_e, \qquad (M_u)
```

Another important additional definition is the `molarmass` constanat, which is automatically selected based on the choice of `boltzmann` constant (but can also be customized if necessary).

```@docs
UnitSystem
```

Other similar packages include [PhysicalConstants.jl](https://github.com/JuliaPhysics/PhysicalConstants.jl), [MathPhysicalConstants.jl](https://github.com/LaGuer/MathPhysicalConstants.jl), [Unitful.jl](https://github.com/PainterQubits/Unitful.jl.git), [UnitfulUS](https://github.com/PainterQubits/UnitfulUS.jl), [UnitfulAstro](https://github.com/JuliaAstro/UnitfulAstro.jl), [UnitfulAtomic](https://github.com/sostock/UnitfulAtomic.jl), [NaturallyUnitful](https://github.com/MasonProtter/NaturallyUnitful.jl), and [UnitfulMoles](https://github.com/rafaqz/UnitfulMoles.jl).

## Metric SI Units

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

Physical measured values with uncertainty are the electron to proton mass ratio `μₑₐ`, proton to atomic mass ratio `μₚₐ`, fine structure constant `αinv`, the Rydberg `R∞` constant, and the Planck mass `mP`.

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

From these numbers along with the `4π*1e-7` value of the Gaussian unit `μ₀`, the constants `planckreduced`, `permeability`, `electronmass`, `molarmass`, and proton to electon mass ratio are computed.

```math
\hbar = \frac{h}{2\pi}, \qquad
\mu_0 = \frac{2h\alpha}{ce^2}, \qquad
m_e = \frac{2hR_\infty}{c\alpha}, \qquad
M_u = \frac{m_e}{\mu_{eu}}N_A = \frac{2h R_\infty N_A}{c\alpha\mu_{eu}}, \qquad
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

Additional reference values include the ground state hyperfine structure transition frequency of caesium-133 `ΔνCs` and luminous efficacy `Kcd` of monochromatic radiation of 540 THz.

```Julia
julia> ΔνCs
9.19263177e9

julia> Kcd
683.0
```

## Other historic systems

Alternatives to the SI unit system are the centimetre-gram-second variants.

```math
\tilde k_B = 10^7\frac{\tilde m_e R_u}{\mu_{eu} \tilde M_u}, \quad
\tilde \hbar = 10^7\hbar, \quad
\tilde c = 100c, \quad
\tilde\mu_0 = 4\pi + \delta\tilde \mu_0, \quad
\tilde m_e = 1000m_e, \quad
(\tilde M_u = 1 + \delta \tilde M_u)
```

```Julia
CGS     ::UnitSystem{1e10*Rᵤ*mₑ/μₑᵤ,1e7*ħ,100𝘤,4π,1000mₑ}
CGS2019 ::UnitSystem{1e7*kB,1e7*ħ,100𝘤,1e7*μ₀,1000mₑ}
```
```@docs
CGS
CGS2019
```

Historically, the `josephson` and `klitzing` constants have been used to define `Conventional` and `CODATA` derived `UnitSystem` variants.

```math
\tilde k_B = \frac{8 \tilde R_\infty R_u}{c\alpha\mu_{eu} \tilde M_u \tilde K_J^2 \tilde R_K}, \quad
\tilde \hbar = \frac{2}{\pi\tilde K_J^2 \tilde R_K}, \quad
\tilde c = c, \quad
\tilde\mu_0 = \frac{2\tilde R_K\alpha}{c}, \quad
\tilde m_e = \frac{8\tilde R_\infty}{\tilde K_J^2\tilde R_Kc\alpha}, \quad
(\tilde M_u = \frac{1}{1000})
```

```Julia
CODATA::UnitSystem{1000Rᵤ*mₑ2014/μₑᵤ,2/RK2014/KJ2014^2/π,𝘤,2RK2014/𝘤/αinv,mₑ2014}()
Conventional::UnitSystem{1000Rᵤ*mₑ/μₑᵤ,2/RK1990/KJ1990^2/π,𝘤,2RK1990/𝘤/αinv,mₑ}()
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

In Britain and the United States an `English` system of engineering units was commonly used.

```math
\tilde k_B = 5.657302466\mathrm{e}{-24}, \quad
\tilde \hbar = \frac{\hbar}{\text{slug}\cdot \text{ft}^2}, \quad
\tilde c = \frac{c}{\text{ft}}, \quad
\tilde\mu_0 = 4\pi, \quad
\tilde m_e = \frac{m_e}{\text{slug}}, \quad
(\tilde M_u = \frac{R_u \tilde m_e ^\circ R}{\tilde k_B\mu_{eu} \text{ft}})
```

```@docs
English
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

## Natural units

With the introduction of the `planckmass` a set of natural atomic unit systems can be derived in terms of the gravitational coupling constant.

```math
\alpha_G = \left(\frac{m_e}{m_P}\right)^2, \qquad
\tilde k_B = 1, \qquad
(\tilde M_u = 1)
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
Schrodinger  ::UnitSystem{1,1,αinv,4π/αinv^2,√αinv*mₑ/mP}
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

## Fundamental constants of physics

The following are fundamental constants of physics:

```math
M_u = m_uN_A = N_A\frac{m_e}{\mu_{eu}} = N_A\frac{m_p}{\mu_{pu}} = N_A\frac{2R_\infty h}{\mu_{eu}c\alpha^2}
```
```@docs
molarmass
```

```math
N_A = \frac{R_u}{k_B} = \frac{M_u}{m_u} = M_u\frac{\mu_{eu}}{m_e} = M_u\frac{\mu_{eu}c\alpha^2}{2R_\infty h}
```
```@docs
avogadro
```

```math
k_B = \frac{R_u}{N_A} = m_u\frac{R_u}{M_u} = \frac{m_e R_u}{\mu_{eu}M_u} = \frac{2R_uR_\infty h}{M_u \mu_{eu}c\alpha^2}
```
```@docs
boltzmann
```

```math
R_u = k_B N_A = \frac{PV}{nT}, \qquad R_s = \frac{R_u}{M_s} = c_p - c_v
```
```@docs
universal
```

```math
c = \frac1{\sqrt{\mu_0\varepsilon_0}}, \qquad \frac{dt}{d\tau} = \sqrt{1-\frac{v^2}{c^2}}
```
```@docs
lightspeed
```

```math
h = 2\pi\hbar = \frac{2e}{K_J} = \frac{8\alpha}{c\mu_0K_J^2} = \frac{4}{K_J^2R_K}
```
```@docs
planck
```

```math
\hbar = \frac{h}{2\pi} = \frac{e}{\pi K_J} = \frac{4\alpha}{\pi c\mu_0K_J^2} = \frac{2}{\pi K_J^2R_K}
```
```@docs
planckreduced
```

```math
m_P = \sqrt{\frac{\hbar c}{G}} = \frac{m_e}{\sqrt{\alpha_G}} = \frac{2R_\infty h}{c\alpha^2\sqrt{\alpha_G}}
```
```@docs
planckmass
```

```math
m_u = \frac{M_u}{N_A} = \frac{m_e}{\mu_{eu}} = \frac{m_p}{\mu_{pu}} = \frac{2R_\infty h}{\mu_{eu}c\alpha^2} = \frac{m_P}{\mu_{eu}}\sqrt{\alpha_G}
```
```@docs
atomicmass
```

```math
m_p = \mu_{pu} m_u = \mu_{pu}\frac{M_u}{N_A} = \mu_{pe}m_e = \mu_{pe}\frac{2R_\infty h}{c\alpha^2} = m_P\mu_{pe}\sqrt{\alpha_G}
```
```@docs
protonmass
```

```math
m_e = \mu_{eu}m_u = \mu_{eu}\frac{M_u}{N_A} = \frac{m_p}{\mu_{pe}} = \frac{2R_\infty h}{c\alpha^2} = m_P\sqrt{\alpha_G}
```
```@docs
electronmass
```

```math
G = \frac{\hbar c}{m_P^2} = \frac{\hbar c\alpha_G}{m_e^2} = \frac{c^3\alpha^2\alpha_G}{8\pi R_\infty^2 h} = \frac{\kappa c^4}{8\pi}
```
```@docs
newton
```

```math
\kappa = \frac{8\pi G}{c^4} = \frac{8\pi\hbar}{c^3m_P^2} = \frac{8\pi\hbar\alpha_G}{c^3m_e^2} = \frac{\alpha^2\alpha_G}{R_\infty^2 h c}
```
```@docs
einstein
```

```math
\sigma = \frac{2\pi^5 k_B^4}{15h^3c^2} = \frac{\pi^2 k_B^4}{60\hbar^3c^2} = \frac{32\pi^5 h}{15c^6\alpha^8} \left(\frac{R_uR_\infty}{\mu_{eu}M_u}\right)^4
```
```@docs
stefan
```

```math
a = 4\frac{\sigma}{c} = \frac{8\pi^5 k_B^4}{15h^3c^3} = \frac{\pi^2 k_B^4}{15\hbar^3c^3} = \frac{2^7\pi^5 h}{15c^7\alpha^8} \left(\frac{R_uR_\infty}{\mu_{eu}M_u}\right)^4
```
```@docs
radiationdensity
```

```math
\mu_0 = \frac{1}{\varepsilon_0 c^2} = \frac{4\pi k_e}{c^2} = \frac{2h\alpha}{ce^2} = \frac{2R_K\alpha}{c}
```
```@docs
permeability
```

```math
\varepsilon_0 = \frac{1}{\mu_0c^2} = \frac{1}{4\pi k_e} = \frac{e^2}{2\alpha hc} = \frac{1}{2R_K\alpha c}
```
```@docs
permittivity
```

```math
k_e = \frac{1}{4\pi\varepsilon_0} = \frac{\mu_0c^2}{4\pi} = \frac{\alpha h c}{2\pi e^2} = \frac{R_K\alpha c}{2\pi}
```
```@docs
coulomb
```

```math
e = \sqrt{\frac{2h\alpha}{Z_0}} = \frac{2}{K_JR_K} = \sqrt{\frac{h}{R_K}} = \frac{hK_J}{2} = \frac{F}{N_A}
```
```@docs
charge
```

```math
F = eN_A = N_A\sqrt{\frac{2h\alpha}{Z_0}} = \frac{2N_A}{K_JR_K} = N_A\sqrt{\frac{h}{R_K}} = \frac{hK_JN_A}{2}
```
```@docs
faraday
```

```math
Z_0 = \mu_0c = \frac{1}{\varepsilon_0 c} = \sqrt{\frac{\mu_0}{\varepsilon_0}} = \frac{2h\alpha}{e^2} = 2R_K\alpha
```
```@docs
impedance
```

```math
G_0 = \frac{2e^2}{h} = \frac{4\alpha}{Z_0} = \frac{2}{R_K} = \frac{hK_J^2}{2} = \frac{2F^2}{hN_A^2}
```
```@docs
conductance
```

```math
R_K = \frac{h}{e^2} = \frac{Z_0}{2\alpha} = \frac{2}{G_0} = \frac{4}{hK_J^2} = h\frac{N_A^2}{F^2}
```
```@docs
klitzing
```

```math
K_J = \frac{2e}{h} = \sqrt{\frac{8\alpha}{hZ_0}} = \sqrt{\frac{4}{hR_K}} = \frac{1}{\Phi_0} = \frac{2F}{hN_A}
```
```@docs
josephson
```

```math
\Phi_0 = \frac{h}{2e} = \sqrt{\frac{hZ_0}{8\alpha}} = \sqrt{\frac{hR_K}{4}} = \frac{1}{K_J} = \frac{hN_A}{2F}
```
```@docs
magneticflux
```

```math
\mu_B = \frac{e\hbar}{2m_e} = \frac{\hbar}{m_eK_JR_K} = \frac{h^2K_J}{8\pi m_e} = \frac{\hbar F}{2m_e N_A} = \frac{ec\alpha^2}{8\pi R_\infty}
```
```@docs
magneton
```

```math
E_h = m_e(c\alpha)^2 = \frac{\hbar c\alpha}{a_0} = \frac{\hbar^2}{m_ea_0^2} = 2R_\infty hc = m_P\sqrt{\alpha_G}(c\alpha)^2
```
```@docs
hartree
```

```math
R_\infty = \frac{E_h}{2hc} = \frac{m_e c\alpha^2}{2h} = \frac{\alpha}{4\pi a_0} = \frac{m_e c r_e}{2ha_0} = \frac{\alpha^2m_ec}{4\pi\hbar}  = \frac{m_Pc\sqrt{\alpha_G}\alpha^2}{2h}
```
```@docs
rydberg
```

```math
a_0 = \frac{\hbar}{m_ec\alpha} = \frac{\hbar^2}{k_e m_ee^2} = \frac{\mu_{pe}}{m_e}a_0^* = \frac{r_e}{\alpha^2} = \frac{\alpha}{4\pi R_\infty}
```
```@docs
bohr
```

```math
a_0* = \frac{m_e}{\mu_{pe}}a_0 = \frac{\hbar}{\mu_{pe}c\alpha} = \frac{\hbar^2}{k_e \mu_{pe}e^2} = \frac{m_er_e}{\mu_{pe}\alpha^2} = \frac{m_e\alpha}{4\pi\mu_{pe} R_\infty}
```
```@docs
bohrreduced
```

```math
r_e = \frac{\hbar\alpha}{m_ec} = \alpha^2a_0 = \frac{e^2 k_e}{m_ec^2} = \frac{2hR_\infty a_0}{m_ec} = \frac{\alpha^3}{4\pi R_\infty}
```
```@docs
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

