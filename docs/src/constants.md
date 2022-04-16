# Physics Constants

```@contents
Pages = ["unitsystems.md"]
Depth = 1
```
```@contents
Pages = ["constants.md"]
```
```@contents
Pages = ["convert.md","units.md"]
Depth = 1
```

The following are fundamental constants of physics:  [![DOI](https://zenodo.org/badge/317419353.svg)](https://zenodo.org/badge/latestdoi/317419353)

```math
\alpha = \frac{\lambda e^2}{4\pi\varepsilon_0\hbar c} = 
\frac{\lambda c\mu_0 (e\alpha_L)^2}{4\pi\hbar} = 
\frac{e^2k_e}{\hbar c} = 
\frac{\lambda e^2}{2\mu_0ch} = 
\frac{\lambda c\mu_0\alpha_L^2}{2R_K} = 
\frac{e^2Z_0}{2h}
```

There exists a deep relationship between the fundamental constants, which also makes them very suitable as a basis for `UnitSystem` dimensional analysis. All of the formulas on this page are part of the `Test` suite to [guarantee](https://github.com/chakravala/UnitSystems.jl/blob/master/test/runtests.jl) their universal correctness.

```math
\mu_{eu} = \frac{m_e}{m_u}, \qquad
\mu_{pu} = \frac{m_p}{m_u}, \qquad
\mu_{pe} = \frac{m_p}{m_e}, \qquad
\alpha_\text{inv} = \frac{1}{\alpha}, \qquad
\alpha_G = \left(\frac{m_e}{m_P}\right)^2
```

```@docs
Universe
```

## Relativistic Constants

```math
c = \frac1{\alpha_L\sqrt{\mu_0\varepsilon_0}} = \frac{1}{\alpha}\sqrt{E_h\frac{g_0}{m_e}} = \frac{g_0\hbar\alpha}{m_e r_e}  = \frac{e^2k_e}{\hbar\alpha} = \frac{m_e^2G}{\hbar\alpha_G}
```
```@docs
MeasureSystems.lightspeed
```

```math
h = 2\pi\hbar = \frac{2e\alpha_L}{K_J} = \frac{8\alpha}{\lambda c\mu_0K_J^2} = \frac{4\alpha_L^2}{K_J^2R_K}
```
```@docs
MeasureSystems.planck
```

```math
\hbar = \frac{h}{2\pi} = \frac{e\alpha_L}{\pi K_J} = \frac{4\alpha}{\pi\lambda c\mu_0K_J^2} = \frac{2\alpha_L}{\pi K_J^2R_K}
```
```@docs
MeasureSystems.planckreduced
```

```math
m_P = \sqrt{\frac{\hbar c}{G}} = \frac{m_e}{\sqrt{\alpha_G}} = \frac{2R_\infty hg_0}{c\alpha^2\sqrt{\alpha_G}}
```
```@docs
planckmass
```

```math
G = k^2 = \frac{\hbar c}{m_P^2} = \frac{\hbar c\alpha_G}{m_e^2} = \frac{c^3\alpha^4\alpha_G}{8\pi g_0^2 R_\infty^2 h} = \frac{\kappa c^4}{8\pi}
```
```@docs
MeasureSystems.gravitation
```

```math
k = \frac{\sqrt{\hbar c}}{m_P} = \frac{\sqrt{\hbar c\alpha_G}}{m_e} = \frac{\alpha^2}{2g_0R_\infty}\sqrt{\frac{c^3\alpha_G}{2\pi h}} = c^2\sqrt{\frac{\kappa}{8\pi}}
```
```@docs
gaussgravitation
```

```math
\kappa = \frac{8\pi G}{c^4} = \frac{8\pi\hbar}{c^3m_P^2} = \frac{8\pi\hbar\alpha_G}{c^3m_e^2} = \frac{\alpha^4\alpha_G}{g_0^2R_\infty^2 h c}
```
```@docs
einstein
```

## Atomic & Nuclear Constants

```math
m_u = \frac{M_u}{N_A} = \frac{m_e}{\mu_{eu}} = \frac{m_p}{\mu_{pu}} = \frac{2R_\infty hg_0}{\mu_{eu}c\alpha^2} = \frac{m_P}{\mu_{eu}}\sqrt{\alpha_G}
```
```@docs
atomicmass
```

```math
m_p = \mu_{pu} m_u = \mu_{pu}\frac{M_u}{N_A} = \mu_{pe}m_e = \mu_{pe}\frac{2R_\infty hg_0}{c\alpha^2} = m_P\mu_{pe}\sqrt{\alpha_G}
```
```@docs
protonmass
```

```math
m_e = \mu_{eu}m_u = \mu_{eu}\frac{M_u}{N_A} = \frac{m_p}{\mu_{pe}} = \frac{2R_\infty h g_0}{c\alpha^2} = m_P\sqrt{\alpha_G}
```
```@docs
MeasureSystems.electronmass
```

```math
E_h = \frac{m_e}{g_0}(c\alpha)^2 = \frac{\hbar c\alpha}{a_0} = \frac{g_0\hbar^2}{m_ea_0^2} = 2R_\infty hc = \frac{m_P}{g_0}\sqrt{\alpha_G}(c\alpha)^2
```
```@docs
hartree
```

```math
R_\infty = \frac{E_h}{2hc} = \frac{m_e c\alpha^2}{2hg_0} = \frac{\alpha}{4\pi a_0} = \frac{m_e r_e c}{2ha_0g_0} = \frac{\alpha^2m_ec}{4\pi\hbar g_0}  = \frac{m_Pc\alpha^2\sqrt{\alpha_G}}{2hg_0}
```
```@docs
rydberg
```

```math
a_0 = \frac{g_0\hbar}{m_ec\alpha} = \frac{g_0\hbar^2}{k_e m_ee^2} = \frac{r_e}{\alpha^2} = \frac{\alpha}{4\pi R_\infty}
```
```@docs
bohr
```

```math
r_e = g_0\frac{\hbar\alpha}{m_ec} = \alpha^2a_0 = g_0\frac{e^2 k_e}{m_ec^2} = \frac{2hR_\infty g_0a_0}{m_ec} = \frac{\alpha^3}{4\pi R_\infty}
```
```@docs
electronradius
```

```math
\Delta\nu_{\text{Cs}} = \Delta\tilde\nu_{\text{Cs}}c = \frac{\Delta\omega_{\text{Cs}}}{2\pi}  = \frac{c}{\Delta\lambda_{\text{Cs}}} = \frac{\Delta E_{\text{Cs}}}{h}
```
```@docs
hyperfine
```

## Thermodynamic Constants

```math
M_u = m_uN_A = N_A\frac{m_e}{\mu_{eu}} = N_A\frac{m_p}{\mu_{pu}} = N_A\frac{2R_\infty hg_0}{\mu_{eu}c\alpha^2}
```
```@docs
MeasureSystems.molarmass
```

```math
N_A = \frac{R_u}{k_B} = \frac{M_u}{m_u} = M_u\frac{\mu_{eu}}{m_e} = M_u\frac{\mu_{eu}c\alpha^2}{2R_\infty h g_0}
```
```@docs
avogadro
```

```math
k_B = \frac{R_u}{N_A} = m_u\frac{R_u}{M_u} = \frac{m_e R_u}{\mu_{eu}M_u} = \frac{2R_uR_\infty h g_0}{M_u \mu_{eu}c\alpha^2}
```
```@docs
MeasureSystems.boltzmann
```

```math
R_u = k_B N_A = k_B\frac{M_u}{m_u} = k_BM_u\frac{\mu_{eu}}{m_e} = k_BM_u\frac{\mu_{eu}c\alpha^2}{2hR_\infty g_0}
```
```@docs
molargas
```

```math
\frac{p_0}{k_B T_0} = \frac{N_Ap_0}{R_uT_0} = \frac{\mu_{eu}M_up_0}{m_e R_u T_0} = \frac{M_u \mu_{eu}c\alpha^2p_0}{2R_uR_\infty hg_0 T_0}
```
```@docs
loschmidt
```

```math
\frac{S_0}{R_u} = log\left(\frac{\hbar^3}{p_0}\sqrt{\left(\frac{m_u}{2\pi g_0}\right)^3 \left(k_BT_0\right)^5}\right)+\frac{5}{2} = log\left(\frac{m_u^4}{p_0}\left(\frac{\hbar}{\sqrt{2\pi g_0}}\right)^3\sqrt{\frac{R_uT_0}{M_u}}^5\right)+\frac{5}{2}
```
```@docs
sackurtetrode
```

```math
\frac{180 R_uV_{it}^2}{43 k_BN_A\Omega_{it}} = 
\frac{180 k_BM_uV_{it}^2}{43 R_um_u\Omega_{it}} = 
\frac{90 k_BM_u\mu_{eu}c\alpha^2V_{it}^2}{43 hg_0R_uR_\infty\Omega_{it}}
```
```@docs
mechanicalheat
```

```math
\sigma = \frac{2\pi^5 k_B^4}{15h^3c^2} = \frac{\pi^2 k_B^4}{60\hbar^3c^2} = \frac{32\pi^5 h}{15c^6\alpha^8} \left(\frac{g_0R_uR_\infty}{\mu_{eu}M_u}\right)^4
```
```@docs
stefan
```

```math
a = 4\frac{\sigma}{c} = \frac{8\pi^5 k_B^4}{15h^3c^3} = \frac{\pi^2 k_B^4}{15\hbar^3c^3} = \frac{2^7\pi^5 h}{15c^7\alpha^8} \left(\frac{g_0R_uR_\infty}{\mu_{eu}M_u}\right)^4
```
```@docs
radiationdensity
```

```math
b = \frac{hc/k_B}{5+W_0(-5 e^{-5})} = \frac{hcM_u/(m_uR_u)}{5+W_0(-5 e^{-5})} = \frac{M_u \mu_{eu}c^2\alpha^2/(2R_uR_\infty g_0)}{5+W_0(-5 e^{-5})}
```
```@docs
wienwavelength
```

```math
\frac{3+W_0(-3 e^{-3})}{h/k_B} = \frac{3+W_0(-3 e^{-3})}{hM_u/(m_uR_u)} = \frac{3+W_0(-3 e^{-3})}{M_u \mu_{eu}c\alpha^2/(2R_uR_\infty g_0)}
```
```@docs
wienfrequency
```

```math
K_{\text{cd}} = \frac{I_v}{\int_0^\infty \bar{y}(\lambda)\cdot\frac{dI_e}{d\lambda}d\lambda}, \qquad 
\bar{y}\left(\frac{c}{540\times 10^{12}}\right)\cdot I_e = 1
```
```@docs
MeasureSystems.luminousefficacy
```

## Electromagnetic Constants

```math
\lambda = \frac{4\pi\alpha_B}{\mu_0\alpha_L} = 4\pi k_e\varepsilon_0 = Z_0\varepsilon_0c
```
```@docs
MeasureSystems.rationalization
```

```math
\alpha_L = \frac{1}{c\sqrt{\mu_0\varepsilon_0}} = \frac{\alpha_B}{\mu_0\varepsilon_0k_e} = \frac{4\pi \alpha_B}{\lambda\mu_0} = \frac{k_m}{\alpha_B}
```
```@docs
MeasureSystems.lorentz
```

```math
\alpha_B = \mu_0\alpha_L\frac{\lambda}{4\pi} = \alpha_L\mu_0\varepsilon_0k_e = \frac{k_m}{\alpha_L} = \frac{k_e}{c}\sqrt{\mu_0\varepsilon_0}
```
```@docs
biotsavart
```

```math
Z_0 = \mu_0\lambda c\alpha_L^2 = \frac{\lambda}{\varepsilon_0 c} = \lambda\alpha_L\sqrt{\frac{\mu_0}{\varepsilon_0}} = \frac{2h\alpha}{e^2} = 2R_K\alpha
```
```@docs
vacuumimpedance
```

```math
\mu_0 = \frac{1}{\varepsilon_0 (c\alpha_L)^2} = \frac{4\pi k_e}{\lambda (c\alpha_L)^2} = \frac{2h\alpha}{\lambda c(e\alpha_L)^2} = \frac{2R_K\alpha}{\lambda c\alpha_L^2}
```
```@docs
MeasureSystems.vacuumpermeability
```

```math
\varepsilon_0 = \frac{1}{\mu_0(c\alpha_L)^2} = \frac{\lambda}{4\pi k_e} = \frac{\lambda e^2}{2\alpha hc} = \frac{\lambda}{2R_K\alpha c}
```
```@docs
vacuumpermittivity
```

```math
k_e = \frac{\lambda}{4\pi\varepsilon_0} = \frac{\mu_0\lambda (c\alpha_L)^2}{4\pi} = \frac{\alpha \hbar c}{e^2} = \frac{R_K\alpha c}{2\pi} = \frac{\alpha_B}{\alpha_L\mu_0\varepsilon_0} = k_mc^2
```
```@docs
electrostatic
```

```math
k_m = \alpha_L\alpha_B = \mu_0\alpha_L^2\frac{\lambda}{4\pi} = \frac{k_e}{c^2} = \frac{\alpha \hbar}{ce^2} = \frac{R_K\alpha}{2\pi c}
```
```@docs
magnetostatic
```

```math
e = \sqrt{\frac{2h\alpha}{Z_0}} = \frac{2\alpha_L}{K_JR_K} = \sqrt{\frac{h}{R_K}} = \frac{hK_J}{2\alpha_L} = \frac{F}{N_A}
```
```@docs
elementarycharge
```

```math
F = eN_A = N_A\sqrt{\frac{2h\alpha}{Z_0}} = \frac{2N_A\alpha_L}{K_JR_K} = N_A\sqrt{\frac{h}{R_K}} = \frac{hK_JN_A}{2\alpha_L}
```
```@docs
faraday
```

```math
G_0 = \frac{2e^2}{h} = \frac{4\alpha}{Z_0} = \frac{2}{R_K} = \frac{hK_J^2}{2\alpha_L^2} = \frac{2F^2}{hN_A^2}
```
```@docs
conductancequantum
```

```math
R_K = \frac{h}{e^2} = \frac{Z_0}{2\alpha} = \frac{2}{G_0} = \frac{4\alpha_L^2}{hK_J^2} = h\frac{N_A^2}{F^2}
```
```@docs
klitzing
```

```math
K_J = \frac{2e\alpha_L}{h} = \alpha_L\sqrt{\frac{8\alpha}{hZ_0}} = \alpha_L\sqrt{\frac{4}{hR_K}} = \frac{1}{\Phi_0} = \frac{2F\alpha_L}{hN_A}
```
```@docs
josephson
```

```math
\Phi_0 = \frac{h}{2e\alpha_L} = \frac{1}{\alpha_L}\sqrt{\frac{hZ_0}{8\alpha}} = \frac{1}{\alpha_L}\sqrt{\frac{hR_K}{4}} = \frac{1}{K_J} = \frac{hN_A}{2F\alpha_L}
```
```@docs
magneticfluxquantum
```

```math
\mu_B = \frac{e\hbar\alpha_L}{2m_e} = \frac{\hbar\alpha_L^2}{m_eK_JR_K} = \frac{h^2K_J}{8\pi m_e} = \frac{\alpha_L\hbar F}{2m_e N_A} = \frac{ec\alpha^2\alpha_L}{8\pi g_0R_\infty}
```
```@docs
magneton
```

## Astronomical Constants

```@docs
solarmass
jupitermass
earthmass
lunarmass
earthradius
greatcircle
gaussianyear
siderealyear
radarmile
hubble
cosmological
```
```math
g_0 = [MLT^{-2}F^{-1}]
```
```@docs
MeasureSystems.gravity
```

## Constants Index

```@index
Pages = ["constants.md","unitsystems.md"]
```

