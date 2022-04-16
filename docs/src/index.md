# Geophysics.jl

*Planetary science data for atmospheric geophysical models*

[![DOI](https://zenodo.org/badge/306497671.svg)](https://zenodo.org/badge/latestdoi/306497671)
[![Build Status](https://travis-ci.org/chakravala/Geophysics.jl.svg?branch=master)](https://travis-ci.org/chakravala/Geophysics.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/dkbkhd26j463hnx7?svg=true)](https://ci.appveyor.com/project/chakravala/geophysics-jl)
[![Coverage Status](https://coveralls.io/repos/chakravala/Geophysics.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/chakravala/Geophysics.jl?branch=master)
[![codecov.io](https://codecov.io/github/chakravala/Geophysics.jl/coverage.svg?branch=master)](https://codecov.io/github/chakravala/Geophysics.jl?branch=master)
[![Gitter](https://badges.gitter.im/Grassmann-jl/community.svg)](https://gitter.im/Grassmann-jl/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)

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

```@contents
Pages = ["unitsystems.md","references.md"]
Depth = 3
```

This package is not limited to atmospheric data: other geophysical data features are intended to be added for oceans, temperature and pressure inside the planets, as well as electrical and magnetic properties of planets.
In this package, any simple Geophysical properties of planets may be added.
Other simple geophysical data about planets, can be added in a collaborative effort.
Complicated models will be excluded from this package, as it is only intended to provide a minimal foundation for geophysical data and constants of various planets, more complicated models should be built separately in packages to build on `Geophysics`.
For example, some geographic conditions can be calculated externally, and then Geophysics is used to load that data.

This `Geophysics` package for the Julia language was created by [github.com/chakravala](https://github.com/chakravala) for research.
