
export Air, AirEnglish, ARDC, US76, ARDCE, US76E
export Standard, Earth1959, Earth1976, Earth1959English, Earth1976English

#   This file is part of Geophysics.jl. It is licensed under the MIT license
#   Geophysics Copyright (C) 2020 Michael Reed

const Air = Gas(0.028965923,720,1008,1.7894e-5,110.4,0.02531,194,288.16) #287.0429-287.058, 0.0289654
const AirEnglish = Gas(0.028965923MEnglish,4290,6006,3.7373e-7,198.72,0.02531kEnglish,349.2,518.69) #1716.49-.5
const MarsAir = Air # to-do

# Atmospheric temperature models

const ARDC = Atmosphere{Air}(
    Values(-6.5e-3,0e-3,3e-3,0e-3,-4.5e-3,0e-3,4e-3),
    Values(-00e3,11e3,25e3,47e3,53e3,79e3,90e3))
const US76 = Atmosphere{Air}(
    Values(-6.5e-3,0.,1e-3,2.8e-3,0e3,-2.8e-3,-2e-3),
    Values(-00e3,11e3,20e3,32e3,47e3,51e3,71e3))
const ARDCE = Atmosphere{AirEnglish}(
    Values(-3.5658e-3,0e-3,1.64584e-3,0e-3,-2.46876e-3,0e-3,-2.19456e-3),
    Values(-0e3,36.089e3,82.021e3,154.199e3,173.885e3,259.176e3,295.276))
const US76E = Atmosphere{AirEnglish}(
    Values(-3.5658e-3,0e3,0.54864e-3,1.53612e-3,0e3,-1.53612e-3,-1.09728e-3),
    Values(-0e3,36.089e3,65.617e3,104.987e3,154.199e3,167.323e3,232.940e3))
const MarsAtmosphere = Atmosphere{MarsAir}(
    Values(-9.88e-4,-2.22e-3),
    Values(-0e3,7e3))

const layers = Values("Troposphere","Tropopause","Stratosphere","Stratosphere","Stratopause","Mesosphere","Mesosphere","Mesopause")

# Standard atmosphere weather conditions (sea level 45° lat)

const Earth1959 = ARDC(288.16)
const Earth1976 = US76(288.15)
const Earth1959English = ARDCE(518.69,2116.2,2.085553e7,32.174)
const Earth1976English = US76E(518.67,2116.2,2.085553e7,32.174)
const MarsWeather = MarsAtmosphere(242.15,699.) # to-do
const english = haskey(ENV,"GEOUNITS") && ENV["GEOUNITS"] == "english"
const Standard = if haskey(ENV,"STDATM") && ENV["STDATM"] == "1976"
    english ? Earth1976English : Earth1976
else
    english ? Earth1959English : Earth1959
end
