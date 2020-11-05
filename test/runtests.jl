using Geophysics, Test

st = Standard(0)
@test gravity() == gravity(0)
@test temperature() == temperature(Standard(0))
@test pressure() == pressure(st)
@test density() ≈ density(st)
@test conductivity() == conductivity(st)
@test elasticity() == elasticity(st)
@test viscosity() == viscosity(st)
@test volume() ≈ volume(st)
@test energy() == energy(st)
@test enthalpy() == enthalpy(st)
@test heatcapacity() ≈ heatcapacity(st)
@test diffusivity() ≈ diffusivity(st)
@test prandtl() ≈ prandtl(st)
@test sonicspeed() == sonicspeed(st)
@test impedance() ≈ impedance(st)
