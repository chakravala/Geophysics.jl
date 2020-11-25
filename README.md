# Geophysics.jl

*Planetary science data for atmospheric geophysical models*

[![DOI](https://zenodo.org/badge/306497671.svg)](https://zenodo.org/badge/latestdoi/306497671)
[![Build Status](https://travis-ci.org/chakravala/Geophysics.jl.svg?branch=master)](https://travis-ci.org/chakravala/Geophysics.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/dkbkhd26j463hnx7?svg=true)](https://ci.appveyor.com/project/chakravala/geophysics-jl)
[![Coverage Status](https://coveralls.io/repos/chakravala/Geophysics.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/chakravala/Geophysics.jl?branch=master)
[![codecov.io](https://codecov.io/github/chakravala/Geophysics.jl/coverage.svg?branch=master)](https://codecov.io/github/chakravala/Geophysics.jl?branch=master)

Provides `Atmosphere` models based on Air Research and Development Command `ARDC` and the United States (1922, 1925, 1956, 1959, 1962, 1966, 1976) Standard Atmosphere `US22,US25,US56,US59,US62,US66,US76` available also in English units `US22E,US25E,US56E,US59E,US62E,US66E,US76E`.
Provided the local absolute sea level and gravitational acceleration, the `Weather` can be initialized based on temperature and pressure.
Presets for the `Standard` atmosphere are provided: `Earth1922`, `Earth1925`, `Earth1956`, `Earth1959`, `Earth1962`, `Earth1966`, `Earth1976`, `Earth1922English`, `Earth1925English`, `Earth1956English`, `Earth1959English`, `Earth1962English`, `Earth1966English`, `Earth1976English`.
By default the 1959 model with metric units is used for `Standard` atmosphere, although a different year can be specified with environment variable `STDATM` and the default unit system can be specified with the `GEOUNITS` environment variable.

```julia
julia> using Geophysics

julia> h = 1000 # altitude, m
1000

julia> gravity(h)
9.803565306802405

julia> temperature(h)
281.66102237169474

julia> pressure(h)
89876.28158431675

julia> sonicspeed(h)
336.4347118683662
```

Values which can be obtained at geometric altitude include `gravity`, `temperature`, `pressure`, `density`, `sonicspeed`, `conductivity`, `viscosity`, `kinematic`, `volume`, `energy`, `enthalpy`, `heatcapacity`, `diffusivity`, `prandtl`, and `impedance`.
In the future, more varieties of atmospheric models will be added for various planets along with winds aloft and turbulent gust distribution data.
Weather data from internet sources may be imported in the future.

This package is not limited to atmospheric data: other geophysical data features are intended to be added for oceans, temperature and pressure inside the planets, as well as electrical and magnetic properties of planets.
In this package, any simple Geophysical properties of planets may be added.
Other simple geophysical data about planets, can be added in a collaborative effort.
Complicated models will be excluded from this package, as it is only intended to provide a minimal foundation for geophysical data and constants of various planets, more complicated models should be built separately in packages to build on `Geophysics`.
For example, some geographic conditions can be calculated externally, and then Geophysics is used to load that data.

## Unit systems

A `UnitSystem` is a consistent set of dimensional values selected to accomodate a particular use case or standardization. In total, five fundamental constants `kB,ħ,𝘤,μ,mₑ` are used to specify a specific unit system. These are the constants of `boltzmann`, `planckreduced`, `lightspeed`, `permeability`, and `electronmass`. Different choices of natural units or physical measurements result in a variety of unit systems optimized for many purposes.

```Julia
    UnitSystem{kB,ħ,𝘤,μ₀,mₑ}
```

Standardized for engineering based on fundamental constants: `kB` Boltzmann's constant, `ħ` reduced Planck's constant, `𝘤` speed of light, `μ₀` vacuum permeability, and `mₑ` electron rest mass.
Primarily the `Metric` SI unit system is used in addition to the historic `English` engineering unit system.
These constants induce derived values for `avogadro`, `boltzmann`, `universal`, `planck`, `planckreduced`, `lightspeed`, `planckmass`, `atomicmass`, `protonmass`, `electronmass`, `newton`, `einstein`, `permeability`, `permittivity`, `coulomb`, and
additional constants `stefan`, `radiationintensity`, `impedance`, `charge`, `magneton`, `conductance`, `faraday`, `magneticflux`, `josephson`, `klitzing`, `hartree`, `rydberg`, `bohr`, `bohrreduced`, and `molarmass`.

https://geophysics.crucialflow.com/dev/units

Additional reference `UnitSystem` variants `CGS`, `CGS2019`, `SI2019`, `CODATA`, `Conventional`; along with several natural atomic units based on the fine structure constant `1/αinv` and the gravitational coupling constant `αG` (`Planck`, `PlanckGauss`, `Stoney`, `Hartree`, `Rydberg`, `Schrodinger`, `Electronic`, `Natural`, `NaturalGauss`, `QCD`, `QCDGauss`, and `QCDoriginal`).

## References
* R. A. Minzer and W. S. Ripley, [The ARDC Model Atmosphere, 1956](https://www.cia.gov/library/readingroom/docs/CIA-RDP81-01043R002600070006-6.pdf), ARDC (1956)
* R. A. Minzer, K. S. W. Champion, and H. L. Pond, [The ARDC Model Atmosphere, 1959](https://apps.dtic.mil/dtic/tr/fulltext/u2/229482.pdf), ARDC (1959)
* NASA, USAF, and USWB, [U.S. Standard Atmosphere, 1962](https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19630003300.pdf), ICAO (1962)
* NOAA, NASA, and USAF, [U.S. Standard Atmosphere, 1976](https://apps.dtic.mil/dtic/tr/fulltext/u2/a035728.pdf), NOAA (1976)
