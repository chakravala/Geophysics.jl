# Unit Conversions

```@contents
Pages = ["unitsystems.md","constants.md"]
Depth = 1
```
```@contents
Pages = ["convert.md"]
```
```@contents
Pages = ["units.md"]
Depth = 1
```

Standardized conversion factors for physics units between `UnitSystem` specifications:  [![DOI](https://zenodo.org/badge/317419353.svg)](https://zenodo.org/badge/latestdoi/317419353)

## Kinematic Units

```@docs
MeasureSystems.A
solidangle
MeasureSystems.T
MeasureSystems.L
area
MeasureSystems.volume
MeasureSystems.wavenumber
angularwavenumber
fuelefficiency
MeasureSystems.frequency
angularfrequency
frequencydrift
MeasureSystems.speed
acceleration
jerk
snap
crackle
pop
volumeflow
```

## Mechanical Units

```@docs
inertia
MeasureSystems.mass
massflow
lineardensity
areadensity
MeasureSystems.density
MeasureSystems.specificweight
MeasureSystems.specificvolume
force
specificforce
gravityforce
MeasureSystems.pressure
compressibility
MeasureSystems.viscosity
diffusivity
rotationalinertia
impulse
momentum
angularmomentum
yank
energy
MeasureSystems.specificenergy
action
fluence
power
powerdensity
MeasureSystems.intensity
spectralflux
soundexposure
impedance(::UnitSystem,::UnitSystem)
MeasureSystems.specificimpedance
admittance
compliance
inertance
```

## Electromagnetic Units

The following unit conversions have been verified for CGS `UnitSystem` variants: [reference](https://www.qsl.net/g4cnn/units/units.htm) [information](https://phys.libretexts.org/Bookshelves/Electricity_and_Magnetism/Book%3A_Electricity_and_Magnetism_(Tatum)/17%3A_Magnetic_Dipole_Moment/17.05%3A_Possible_Alternative_Definitions_of_Magnetic_Moment).

```@docs
charge(::UnitSystem,::UnitSystem)
chargedensity
linearchargedensity
exposure
mobility
current
currentdensity
resistance
conductance(::UnitSystem,::UnitSystem)
resistivity
conductivity
capacitance
inductance
reluctance
permeance
permittivity(::UnitSystem,::UnitSystem)
MeasureSystems.permeability(::UnitSystem,::UnitSystem)
susceptibility
specificsusceptibility
demagnetizingfactor
vectorpotential
electricpotential
magneticpotential
electricfield
magneticfield
electricflux
magneticflux(::UnitSystem,::UnitSystem)
electricfluxdensity
magneticfluxdensity
electricdipolemoment
magneticdipolemoment
electricpolarizability
magneticpolarizability
magneticmoment
specificmagnetization
polestrength
```

## Thermodynamic Units

```@docs
MeasureSystems.temperature
entropy
specificentropy
volumeheatcapacity
MeasureSystems.thermalconductivity
thermalconductance
thermalresistivity
thermalresistance
thermalexpansion
lapserate
```

## Molar Units

```@docs
MeasureSystems.molarmass(::UnitSystem,::UnitSystem)
molality
molaramount
molarity
molarvolume
molarentropy
molarenergy
molarconductivity
molarsusceptibility
catalysis
specificity
```

## Photometric Units

```@docs
luminousflux
luminousintensity
luminance
illuminance
luminousenergy
luminousexposure
MeasureSystems.luminousefficacy(::UnitSystem,::UnitSystem)
```

## Conversion Index

```@index
Pages = ["convert.md","unitsystems.md"]
```
