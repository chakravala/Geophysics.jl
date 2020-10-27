# Geophysics.jl

*Planetary science data for atmospheric geophysical models*



[![DOI](https://zenodo.org/badge/306497671.svg)](https://zenodo.org/badge/latestdoi/306497671)
[![Build Status](https://travis-ci.org/chakravala/Geophysics.jl.svg?branch=master)](https://travis-ci.org/chakravala/Geophysics.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/dkbkhd26j463hnx7?svg=true)](https://ci.appveyor.com/project/chakravala/geophysics-jl)
[![Coverage Status](https://coveralls.io/repos/chakravala/Geophysics.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/chakravala/Geophysics.jl?branch=master)
[![codecov.io](https://codecov.io/github/chakravala/Geophysics.jl/coverage.svg?branch=master)](https://codecov.io/github/chakravala/Geophysics.jl?branch=master)

Provides `Atmosphere` models based on Air Research and Development Command `ARDC` and the United States 1976 Standard Atmosphere `US76` available also in English units `ARDCE` and `US76E`.
Provided the local absolute sea level and gravitational acceleration, the `Weather` can be initialized based on temperature and pressure.
Presets for the `Standard` atmosphere are provided: `Earth1959`, `Earth1976`, `Earth1959English`, `Earth1976English`.
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

Values which can be obtained at geometric altitude include `gravity`, `temperature`, `pressure`, `density`, and `sonicspeed`.
In the future, more varieties of atmospheric models will be added for various planets along with winds aloft and turbulent gust distribution data.
Weather data from internet sources may be imported in the future.
Additionally this package is not limited to atmospheric data: other geophysical data features can be added for oceans, temperature and pressure inside the planets, as well as electrical and magnetic properties of planets.
In this package, any simple Geophysical properties of planets may be added.
Other simple geophysical data about planets, can be added in a collaborative effort.
Complicated models will be excluded from this package, as it is only intended to provide a minimal foundation for geophysical data and constants of various planets, more complicated models should be built separately in packages to build on `Geophysics`.
For example, some geographic conditions can be calculated externally, and then Geophysics is used to load that data.

## References
* R. A. Minzer, K. S. W. Champion, and H. L. Pond, [The ARDC Model Atmosphere](https://apps.dtic.mil/dtic/tr/fulltext/u2/229482.pdf), ARDC (1959)
* NOAA, NASA, and USAF, [U.S. Standard Atmosphere 1976](https://apps.dtic.mil/dtic/tr/fulltext/u2/a035728.pdf), NOAA (1976)
