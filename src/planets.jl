
export Air, AirEnglish, ARDC, ARDCE, Standard
export US22, US25, US56, US59, US62, US66, US76
export US22E, US25E, US56E, US59E, US62E, US66E, US76E
export Earth1922, Earth1925, Earth1956, Earth1959, Earth1962, Earth1966, Earth1976
export Earth1922English, Earth1925English, Earth1956English
export Earth1959English, Earth1962English, Earth1966English, Earth1976English

#   This file is part of Geophysics.jl. It is licensed under the AGPL license
#   Geophysics Copyright (C) 2020 Michael Reed

export Nitrogen, Oxygen, CarbonDioxide, Methane, Hydrogen
export Argon, Neon, Helium, Krypton, Xenon
export N2, O2, CO2, CH4, H2, N₂, O₂, CO₂, CH₄, H₂, Ar, Ne, He, Kr, Xe

export NitrogenEnglish, OxygenEnglish, CarbonDioxideEnglish, MethaneEnglish, HydrogenEnglish
export ArgonEnglish, NeonEnglish, HeliumEnglish, KryptonEnglish, XenonEnglish
export N2E, O2E, CO2E, CH4E, H2E, N₂E, O₂E, CO₂E, CH₄E, H₂E, ArE, NeE, HeE, KrE, XeE

export AirMix, Nitrox, Traces
export AirMixEnglish, NitroxEnglish, TracesEnglish

const Nitrogen = DiatomicGas(28.013e-3,2744e2,1.735e-5,107,25.11e-3,150)
const Oxygen = DiatomicGas(31.999e-3,2061e2,1.999e-5,139,25.33e-3,240)
const Argon = AtomicGas(39.948e-3,2.187e-5,144,17.23e-3,170)
const CarbonDioxide = TriatomicGas(44.01e-3,2565e2,1480e2,14.45e-5,222,15.8e-3,1800)
const Neon = AtomicGas(20.18e-3,3.078e-5,NaN,48.29e-3,NaN)
const Helium = AtomicGas(4.003e-3,2.928e-5,NaN,152.07e-3,NaN)
const Methane = PentatomicGas(16.042e-3,10.74e-5,NaN,32.7e-3,NaN)
const Krypton = AtomicGas(82.798e-3,2.432e-5,NaN,9.12e-3,NaN)
const Hydrogen  = DiatomicGas(2.016e-3,4342e2,0.866e-5,97,180.1e-3,120)
const Xenon = AtomicGas(131.293e-3,2.229e-5,NaN,5.27e-3,NaN)
const N2,O2,CO2,CH4,H2 = Nitrogen,Oxygen,CarbonDioxide,Methane,Hydrogen
const N₂,O₂,CO₂,CH₄,H₂ = Nitrogen,Oxygen,CarbonDioxide,Methane,Hydrogen
const Ar,Ne,He,Kr,Xe = Argon,Neon,Helium,Krypton,Xenon

const NitrogenEnglish = DiatomicGas(Nitrogen)
const OxygenEnglish = DiatomicGas(Oxygen)
const ArgonEnglish = AtomicGas(Argon)
const CarbonDioxideEnglish = TriatomicGas(CarbonDioxide)
const NeonEnglish = AtomicGas(Neon)
const HeliumEnglish = AtomicGas(Helium)
const MethaneEnglish = PentatomicGas(Methane)
const KryptonEnglish = AtomicGas(Krypton)
const HydrogenEnglish = DiatomicGas(Hydrogen)
const XenonEnglish = AtomicGas(Xenon)
const N2E,O2E,CO2E,CH4E,H2E = NitrogenEnglish,OxygenEnglish,CarbonDioxideEnglish,MethaneEnglish,HydrogenEnglish
const N₂E,O₂E,CO₂E,CH₄E,H₂E = NitrogenEnglish,OxygenEnglish,CarbonDioxideEnglish,MethaneEnglish,HydrogenEnglish
const ArE,NeE,HeE,KrE,XeE = ArgonEnglish,NeonEnglish,HeliumEnglish,KryptonEnglish,XenonEnglish

#*Nitrogen N2: 78.084%, 28.013, 1.735, 25.11
#*Oxygen O2: 20.946%, 31.999, 1.999, 25.33
#*Argon Ar: 0.934%, 39.948, 2.187, 17.23
#*CO2: 0.033%, 44.01, 14.45, 15.8
# Neon Ne: 0.001818%, 20.18, 3.078, 48.29
# Helium He: 0.000524%, 4.003, 2.928, 152.07
# Methane CH4: 0.000179%, 16.042, 10.74, 32.7
# Krypton Kr: 0.0001%, 82.798, 2.432, 9.12
#*Hydrogen H2: 0.00005%, 2.016, 0.866, 180.1
# Xenon Xe: 0.000009%, 131.293, 2.229, 5.27
# Water H20: 0%,

const air = SutherlandGas(0.028965923,720,1.7894e-5,110.4,0.02531,194,288.16) #287.0429-287.058, 0.0289654
const airEnglish = SutherlandGas(air,4290) # 3.7373e-7, 1716.49-.5
#const MarsAir = Air # to-do
export air, airEnglish

const Air,AirEnglish = if !(VERSION<v"1.3")
#const Nitrox = 0.761N₂+0.239O₂ #+0.0000018078CH₄
const Nitrox = 0.7808093N₂+0.2094552O₂+0.009338Ar+0.0003975CO₂
const Traces = 0.726803Ne+0.20966He+0.04006Kr+0.019824H₂+0.003653Xe
const MainGases = 0.7808089N₂+0.2094551O₂+0.009338Ar+0.0003976CO₂+4.96e-7*H₂
const TraceGases = 0.7415Ne+0.2139He+0.040871Kr+0.003729Xe
const AirMix = 0.7807898N₂+0.20945O₂+0.0093378Ar+0.00039738CO₂+0.000018186Ne+0.0000052461He+0.0000010024Kr+0.00000049603H₂+0.000000091399Xe
#const NitroxEnglish = 0.761N₂E+0.239O₂E #+0.0000018078CH₄E
const NitroxEnglish = 0.7808093N₂E+0.2094552O₂E+0.009338ArE+0.0003975CO₂E
const TracesEnglish = 0.726803NeE+0.20966HeE+0.04006KrE+0.019824H₂E+0.003653XeE
const AirMixEnglish = 0.7807898N₂E+0.20945O₂E+0.0093378ArE+0.00039738CO₂E+0.000018186NeE+0.0000052461HeE+0.0000010024KrE+0.00000049603H₂E+0.000000091399XeE
Nitrox,NitroxEnglish else air,airEnglish end

# Atmospheric temperature models

const US22 = Atmosphere{Air}(Values(-6.5e-3,0e-3),Values(-00e3,11e3))
const US25 = Atmosphere{Air}(Values(-6.5e-3,0e-3),Values(-00e3,10.76923e3))
const US56 = Atmosphere{Air}(
    Values(-6.5e-3,0e-3,3e-3,0e-3,-3.9e-3,0e-3,3.5e-3,10e-3,5.8e-3),
    Values(-00e3,11e3,25e3,47e3,53e3,75e3,90e3,126e3,175e3))
const US59 = Atmosphere{Air}(
    Values(-6.5e-3,0e-3,3e-3,0e-3,-4.5e-3,0e-3,4e-3,20e-3,10e-3,5e-3,3.5e-3),
    Values(-00e3,11e3,25e3,47e3,53e3,79e3,90e3,105e3,160e3,170e3,200e3))
#=const US62 = Atmosphere{Air}(
    Values(-6.5e-3,0.,1e-3,2.8e-3,0e3,-2e-3,-4e-3,0.,3e-3),
    Values(-00e3,11e3,20e3,32e3,47e3,52e3,61e3,79e3,90e3))=#
const US62 = Atmosphere{Air}(
    Values(-6.5e-3,0.,1e-3,2.8e-3,0e3,-2e-3,-4e-3,0.,3e-3, 5e-3, 10e-3, 20e-3, 15e-3, 10e-3, 7e-3, 5e-3, 4e-3, 3.3e-3, 2.6e-3, 1.7e-3, 1.1e-3),
    Values(-00e3,11e3,20e3,32e3,47e3,52e3,61e3,79e3,90e3, 100e3, 110e3, 120e3, 150e3, 160e3, 170e3, 190e3, 230e3, 300e3, 400e3, 600e3, 700e3))
const US66 = Atmosphere{Air}(
    Values(-6.5e-3,0.,1e-3,2.8e-3,0e3,-2e-3,-3.9e-3,0.,3e-3),
    Values(-00e3,11e3,20.1e3,32.2e3,47.3e3,52.4e3,61.6e3,80e3,90e3))
const US76 = Atmosphere{Air}(
    Values(-6.5e-3,0.,1e-3,2.8e-3,0e3,-2.8e-3,-2e-3,0.,-Inf,12e-3,Inf),
    Values(-00e3,11e3,20e3,32e3,47e3,51e3,71e3,86e3,91e3,110e3,120e3))
const US22E = Atmosphere{AirEnglish}(Values(-3.5658e-3,0e-3),Values(-0e3,36.089e3))
const US25E = Atmosphere{AirEnglish}(Values(-3.5658e-3,0e-3),Values(-0e3,35.332e3))
const US56E = Atmosphere{AirEnglish}(
    Values(-3.5658e-3,0.,1.64584e-3,0.,-2.1397e-3,0.,1.92024e-3,5.4864e-3,3.1821e-3),
    Values(-0.,36089.,82021.,154199.,173885.,246063.,295276.,413386.,574147.))
const US59E = Atmosphere{AirEnglish}(
    Values(-3.5658e-3,0.,1.64584e-3,0.,-2.46876e-3,0.,2.19456e-3,10.9728e-3,5.4864e-3,2.7432e-3,1.92024e-3),
    Values(-0.,36089.,82021.,154199.,173885.,259176.,295276.,344488.,524934.,557743.,656168.))
const US62E = Atmosphere{AirEnglish}(
    Values(-3.5658e-3,0.,0.54864e-3,1.53612e-3,0e3,-1.09728e-3,-2.1946e-3,0.,1.6459e-3),
    Values(-0.,36089.,65617.,104987.,154199.,170604.,200131.,259186.,295276.))
const US66E = Atmosphere{AirEnglish}(
    Values(-3.5658e-3,0.,0.54864e-3,1.53612e-3,0e3,-1.09728e-3,-2.1397e-3,0.,1.6459e-3),
    Values(-0.,36089.,65945.,105643.,155184.,171916.,202.1e3,262467.,295276.))
const US76E = Atmosphere{AirEnglish}(
    Values(-3.5658e-3,0e3,0.54864e-3,1.53612e-3,0e3,-1.53612e-3,-1.09728e-3),
    Values(-0e3,36.089e3,65.617e3,104.987e3,154.199e3,167.323e3,232.940e3))
#=const MarsAtmosphere = Atmosphere{MarsAir}(
    Values(-9.88e-4,-2.22e-3),
    Values(-0e3,7e3))=#

const ARDC,ARDCE = US59,US59E # deprecate?

const layers = Values("Troposphere","Tropopause","Stratosphere","Stratosphere","Stratopause","Mesosphere","Mesosphere","Mesopause")

# Standard atmosphere weather conditions (sea level 45° lat)

const Earth1922,Earth1922English = US22(288.16),US22E(518.69,2116.2,2.085553e7,32.174)
const Earth1925,Earth1925English = US25(288.16),US25E(518.69,2116.2,2.085553e7,32.174)
const Earth1956,Earth1956English = US56(288.16),US56E(518.69,2116.2,2.085553e7,32.174)
const Earth1959,Earth1959English = US59(288.16),US59E(518.69,2116.2,2.085553e7,32.174)
const Earth1962,Earth1962English = US62(288.15),US62E(518.67,2116.2,2.085553e7,32.174)
const Earth1966,Earth1966English = US66(288.15),US66E(518.67,2116.2,2.085553e7,32.174)
const Earth1976,Earth1976English = US76(288.15),US76E(518.67,2116.2,2.085553e7,32.174)
#const MarsWeather = MarsAtmosphere(242.15,699.) # to-do
const english = haskey(ENV,"GEOUNITS") && ENV["GEOUNITS"] == "english"
const Standard = if haskey(ENV,"STDATM")
    if ENV["STDATM"] == "1922"
        english ? Earth1922English : Earth1922
    elseif ENV["STDATM"] == "1925"
        english ? Earth1925English : Earth1925
    elseif ENV["STDATM"] == "1956"
        english ? Earth1956English : Earth1956
    elseif ENV["STDATM"] == "1959"
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
