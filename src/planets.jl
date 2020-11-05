
export Air, AirEnglish, ARDC, US62, US66, US76, ARDCE, US62E, US66E, US76E
export Standard, Earth1959, Earth1976, Earth1959English, Earth1976English
export Earth1962, Earth1966, Earth1962English, Earth1966English

#   This file is part of Geophysics.jl. It is licensed under the MIT license
#   Geophysics Copyright (C) 2020 Michael Reed

const Air = Gas(0.028965923,720,1008,1.7894e-5,110.4,0.02531,194,288.16) #287.0429-287.058, 0.0289654
const AirEnglish = Gas(0.028965923MEnglish,4290,6006,3.7373e-7,198.72,0.02531kEnglish,349.2,518.69) #1716.49-.5
const MarsAir = Air # to-do

# Atmospheric temperature models

const ARDC = Atmosphere{Air}(
    Values(-6.5e-3,0e-3,3e-3,0e-3,-4.5e-3,0e-3,4e-3),
    Values(-00e3,11e3,25e3,47e3,53e3,79e3,90e3))
const US62 = Atmosphere{Air}(
    Values(-6.5e-3,0.,1e-3,2.8e-3,0e3,-2e-3,-4e-3,0.,3e-3),
    Values(-00e3,11e3,20e3,32e3,47e3,52e3,61e3,79e3,90e3))
const US66 = Atmosphere{Air}(
    Values(-6.5e-3,0.,1e-3,2.8e-3,0e3,-2e-3,-3.9e-3,0.,3e-3),
    Values(-00e3,11e3,20.1e3,32.2e3,47.3e3,52.4e3,61.6e3,80e3,90e3))
const US76 = Atmosphere{Air}(
    Values(-6.5e-3,0.,1e-3,2.8e-3,0e3,-2.8e-3,-2e-3,0.),
    Values(-00e3,11e3,20e3,32e3,47e3,51e3,71e3,86e3))
const ARDCE = Atmosphere{AirEnglish}(
    Values(-3.5658e-3,0e-3,1.64584e-3,0e-3,-2.46876e-3,0e-3,-2.19456e-3),
    Values(-0e3,36.089e3,82.021e3,154.199e3,173.885e3,259.176e3,295.276e3))
const US62E = Atmosphere{AirEnglish}(
    Values(-3.5658e-3,0.,0.54864e-3,1.53612e-3,0e3,-1.09728e-3,-2.1946e-3,0.,1.6459e-3),
    Values(-0.,36089.,65617.,104987.,154199.,170604.,200131.,259186.,295276.))
const US66E = Atmosphere{AirEnglish}(
    Values(-3.5658e-3,0.,0.54864e-3,1.53612e-3,0e3,-1.09728e-3,-2.1397e-3,0.,1.6459e-3),
    Values(-0.,36089.,65945.,105643.,155184.,171916.,202.1e3,262467.,295276.))
const US76E = Atmosphere{AirEnglish}(
    Values(-3.5658e-3,0e3,0.54864e-3,1.53612e-3,0e3,-1.53612e-3,-1.09728e-3),
    Values(-0e3,36.089e3,65.617e3,104.987e3,154.199e3,167.323e3,232.940e3))
const MarsAtmosphere = Atmosphere{MarsAir}(
    Values(-9.88e-4,-2.22e-3),
    Values(-0e3,7e3))

const layers = Values("Troposphere","Tropopause","Stratosphere","Stratosphere","Stratopause","Mesosphere","Mesosphere","Mesopause")

# Standard atmosphere weather conditions (sea level 45° lat)

const Earth1959,Earth1959English = ARDC(288.16),ARDCE(518.69,2116.2,2.085553e7,32.174)
const Earth1962,Earth1962English = US62(288.15),US62E(518.67,2116.2,2.085553e7,32.174)
const Earth1966,Earth1966English = US66(288.15),US66E(518.67,2116.2,2.085553e7,32.174)
const Earth1976,Earth1976English = US76(288.15),US76E(518.67,2116.2,2.085553e7,32.174)
const MarsWeather = MarsAtmosphere(242.15,699.) # to-do
const english = haskey(ENV,"GEOUNITS") && ENV["GEOUNITS"] == "english"
const Standard = if haskey(ENV,"STDATM")
    if ENV["STDATM"] == "1959"
        english ? Earth1959English : Earth1959
    elseif ENV["STDATM"] == "1962"
        english ? Earth1962English : Earth1962
    elseif ENV["STDATM"] == "1966"
        english ? Earth1966English : Earth1966
    elseif ENV["STDATM"] == "1976"
        english ? Earth1976English : Earth1976
    else
        throw(error("unsupported STDATM environment"))
    end
else
    english ? Earth1959English : Earth1959
end
