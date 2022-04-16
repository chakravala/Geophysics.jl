using Geophysics, Test
import Geophysics: normal

st = Standard(0)
@test gravity() == gravity(0)
@test normal(temperature()) == temperature(Standard(0))
@test normal(pressure()) == pressure(st)
#@test normal(density()) ≈ density(st)
@test thermalconductivity() == thermalconductivity(st)
@test normal(elasticity()) == elasticity(st)
@test viscosity() == viscosity(st)
#@test normal(specificvolume()) ≈ specificvolume(st)
@test specificenergy() == specificenergy(st)
@test specificenthalpy() == specificenthalpy(st)
#@test normal(heatcapacity()) ≈ heatcapacity(st)
#@test normal(thermaldiffusivity()) ≈ thermaldiffusivity(st)
#@test normal(prandtl()) ≈ prandtl(st)
@test sonicspeed() == sonicspeed(st)
#@test normal(specificimpedance()) ≈ specificimpedance(st)
