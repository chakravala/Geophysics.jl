# Unit Conversions

```@contents
Pages = ["units.md","constants.md","convert.md"]
```

Common conversion factors for physics units between `UnitSystem` specifications.

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
gforce
stiffness
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

**Warning**: the following unit conversions have not yet been verified for CGS `UnitSystem` variants due to lack of [reference](https://www.qsl.net/g4cnn/units/units.htm) [information](https://phys.libretexts.org/Bookshelves/Electricity_and_Magnetism/Book%3A_Electricity_and_Magnetism_(Tatum)/17%3A_Magnetic_Dipole_Moment/17.05%3A_Possible_Alternative_Definitions_of_Magnetic_Moment): `rigidity`, `mobility`, `magneticmoment`.

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
magnetizability
magnetization
specificmagnetization
rigidity
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
thermalresistance
thermalexpansion
lapserate
```

## Molar Units

```@docs
MeasureSystems.molarmass(::UnitSystem,::UnitSystem)
molality
mole
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
luminance
luminousenergy
luminousexposure
MeasureSystems.luminousefficacy(::UnitSystem,::UnitSystem)
```

## Conversion Index

```@index
Pages = ["convert.md","units.md"]
```
