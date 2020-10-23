using Geophysics, Test

@test gravity() == gravity(0)
@test temperature() == temperature(0)
@test pressure() == pressure(0)
@test density() == density(0)
@test sonicspeed() == sonicspeed(0)
