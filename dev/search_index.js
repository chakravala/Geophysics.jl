var documenterSearchIndex = {"docs":
[{"location":"references.html#References-1","page":"References","title":"References","text":"","category":"section"},{"location":"references.html#","page":"References","title":"References","text":"(Image: DOI) (Image: Build Status) (Image: Build status) (Image: Docs Stable) (Image: Docs Dev) (Image: Gitter)","category":"page"},{"location":"references.html#","page":"References","title":"References","text":"R. A. Minzer and W. S. Ripley, The ARDC Model Atmosphere, 1956, ARDC (1956)\nR. A. Minzer, K. S. W. Champion, and H. L. Pond, The ARDC Model Atmosphere, 1959, ARDC (1959)\nNASA, USAF, and USWB, U.S. Standard Atmosphere, 1962, ICAO (1962)\nNOAA, NASA, and USAF, U.S. Standard Atmosphere, 1976, NOAA (1976)","category":"page"},{"location":"units.html#Dimensional-unit-systems-1","page":"Unit systems","title":"Dimensional unit systems","text":"","category":"section"},{"location":"units.html#","page":"Unit systems","title":"Unit systems","text":"Pages = [\"units.md\",\"references.md\"]","category":"page"},{"location":"units.html#","page":"Unit systems","title":"Unit systems","text":"A UnitSystem is a consistent set of dimensional values selected to accomodate a particular use case or standardization. In total, five fundamental constants kB,ħ,𝘤,μ,mₑ are used to specify a specific unit system. These are the constants of boltzmann, planckreduced, lightspeed, permeability, and electronmass. Different choices of natural units or physical measurements result in a variety of unit systems optimized for many purposes.","category":"page"},{"location":"units.html#","page":"Unit systems","title":"Unit systems","text":"UnitSystem","category":"page"},{"location":"units.html#Geophysics.UnitSystem","page":"Unit systems","title":"Geophysics.UnitSystem","text":"UnitSystem{kB,ħ,𝘤,μ₀,mₑ}\n\nStandardized for engineering based on fundamental constants: kB Boltzmann's constant, ħ reduced Planck's constant, 𝘤 speed of light, μ₀ vacuum permeability, and mₑ electron rest mass. Primarily the Metric SI unit system is used in addition to the historic English engineering unit system. These constants induce derived values for avogadro, boltzmann, universal, planck, planckreduced, lightspeed, planckmass, atomicmass, protonmass, electronmass, newton, einstein, permeability, permittivity, coulomb, and additional constants stefan, radiationintensity, impedance, charge, magneton, conductance, faraday, josephson, magneticflux, klitzing, hardtree, rydberg, bohr, and bohrreduced.\n\nAdditional reference UnitSystem variants CGS, CGS2019, SI2019, CODATA, Conventional; along with several natural atomic units based on the fine structure constant 1/αinv and the gravitational coupling constant αG (Planck, PlanckGauss, Stoney, Hartree, Rydberg, Schrodinger, Electronic, Natural, NaturalGauss, QCD, QCDGauss, and QCDoriginal).\n\n\n\n\n\n","category":"type"},{"location":"units.html#Metric-SI-Units-1","page":"Unit systems","title":"Metric SI Units","text":"","category":"section"},{"location":"units.html#","page":"Unit systems","title":"Unit systems","text":"In the Systeme International d'Unites (the SI units) the UnitSystem constants are derived from the most accurate possible physical measurements and a few exactly defined constants. Exact values are the avogadro number, boltzmann constant, planck constant, lightspeed definition, and elementary charge definition.","category":"page"},{"location":"units.html#","page":"Unit systems","title":"Unit systems","text":"julia> NA # avogadro\n6.02214076e23\n\njulia> kB # boltzmann\n1.380649e-23\n\njulia> 𝘩 # planck\n6.62607015e-34\n\njulia> 𝘤 # lightspeed\n2.99792458e8\n\njulia> 𝘦 # charge\n1.602176634e-19","category":"page"},{"location":"units.html#","page":"Unit systems","title":"Unit systems","text":"Physical measured values with uncertainty are the electron to proton mass ratio μₑₐ, proton to atomic mass ratio μₚₐ, fine structure constant αinv, the Rydberg R∞ constant, and the Planck mass mP.","category":"page"},{"location":"units.html#","page":"Unit systems","title":"Unit systems","text":"julia> μₑₐ\n0.0005485799090649074\n\njulia> μₚₐ\n1.007276466621\n\njulia> αinv\n137.035999084\n\njulia> R∞ # rydberg\n1.09737315681601e7\n\njulia> mP # planckmass\n2.176434e-8","category":"page"},{"location":"units.html#","page":"Unit systems","title":"Unit systems","text":"From these numbers along with the 4π*1e-7 value of the Gaussian unit μ₀, the molar mass constant Mᵤ, planckreduced, permeability, electronmass, and proton to electon mass ratio are computed.","category":"page"},{"location":"units.html#","page":"Unit systems","title":"Unit systems","text":"julia> Mᵤ # αinv^2*R∞*NA*2𝘩/𝘤/μₑₐ\n0.000999999999656256\n\njulia> ħ # 𝘩/2π\n1.0545718176461565e-34\n\njulia> μ₀+δμ₀ # 2𝘩/𝘤/αinv/𝘦^2\n1.256637062121048e-6\n\njulia> mₑ # Mᵤ*μₑₐ/NA\n9.109383701558256e-31\n\njulia> μₚₑ # μₚₐ/μₑₐ\n1836.152673432705","category":"page"},{"location":"units.html#","page":"Unit systems","title":"Unit systems","text":"This results in the Metric::UnitSystem{kB,ħ,𝘤,μ₀,mₑ} and SI2019::UnitSystem{kB,ħ,𝘤,μ₀+δμ₀,mₑ} variants.","category":"page"},{"location":"units.html#","page":"Unit systems","title":"Unit systems","text":"Geophysics.Metric\nGeophysics.SI2019","category":"page"},{"location":"units.html#Geophysics.Metric","page":"Unit systems","title":"Geophysics.Metric","text":"Metric\n\nSysteme International d'Unites (the SI units) adopted as the preffered UnitSystem.\n\njulia> boltzmann(Metric) # J⋅K⁻¹\n1.380649e-23\n\njulia> planckreduced(Metric) # J⋅s⋅rad⁻¹\n1.0545718176461565e-34\n\njulia> lightspeed(Metric) # m⋅s⁻¹\n2.99792458e8\n\njulia> permeability(Metric) # H⋅m⁻¹\n1.2566370614359173e-6\n\njulia> electronmass(Metric) # kg\n9.109383701558256e-31\n\n\n\n\n\n","category":"constant"},{"location":"units.html#Geophysics.SI2019","page":"Unit systems","title":"Geophysics.SI2019","text":"SI2019\n\nSysteme International d'Unites (the SI units) with μ₀+6.851306461996397e-16 for a tuned charge.\n\njulia> boltzmann(SI2019) # J⋅K⁻¹\n1.380649e-23\n\njulia> planckreduced(SI2019) # J⋅s⋅rad⁻¹\n1.0545718176461565e-34\n\njulia> lightspeed(SI2019) # m⋅s⁻¹\n2.99792458e8\n\njulia> permeability(SI2019) # H⋅m⁻¹\n1.2566370619358342e-6\n\njulia> electronmass(SI2019) # kg\n9.109383701558256e-31\n\n\n\n\n\n","category":"constant"},{"location":"units.html#","page":"Unit systems","title":"Unit systems","text":"Additional reference values include the ground state hyperfine structure transition frequency of caesium-133 ΔνCs and luminous efficacy Kcd of monochromatic radiation of 540 THz.","category":"page"},{"location":"units.html#","page":"Unit systems","title":"Unit systems","text":"julia> ΔνCs\n9.19263177e9\n\njulia> Kcd\n683.0","category":"page"},{"location":"units.html#Other-historic-systems-1","page":"Unit systems","title":"Other historic systems","text":"","category":"section"},{"location":"units.html#","page":"Unit systems","title":"Unit systems","text":"Alternatives to the SI unit system are the centimetre-gram-second variants.","category":"page"},{"location":"units.html#","page":"Unit systems","title":"Unit systems","text":"CGS     ::UnitSystem{1e7*kB,1e7*ħ,100𝘤,1e7*μ₀,1000mₑ}\nCGS2019 ::UnitSystem{1e7*kB,1e7*ħ,100𝘤,1e7(μ₀+δμ₀),1000mₑ}","category":"page"},{"location":"units.html#","page":"Unit systems","title":"Unit systems","text":"Geophysics.CGS\nGeophysics.CGS2019","category":"page"},{"location":"units.html#Geophysics.CGS","page":"Unit systems","title":"Geophysics.CGS","text":"CGS\n\nCentimetre-gram-second UnitSystem variant of Metric system based on factors of 1e2,1e3.\n\njulia> boltzmann(CGS) # erg⋅K⁻¹\n1.380649e-16\n\njulia> planckreduced(CGS) # erg⋅s⋅rad⁻¹\n1.0545718176461565e-27\n\njulia> lightspeed(CGS) # cm⋅s⁻¹\n2.99792458e8\n\njulia> permeability(CGS) # erg⋅A⁻²⋅cm⁻¹\n12.566370614359172\n\njulia> electronmass(CGS) # g\n9.109383701558256e-28\n\n\n\n\n\n","category":"constant"},{"location":"units.html#Geophysics.CGS2019","page":"Unit systems","title":"Geophysics.CGS2019","text":"CGS2019\n\nCentimetre-gram-second UnitSystem variant of the tuned SI2019 unit specification.\n\njulia> boltzmann(CGS2019) # erg⋅K⁻¹\n1.380649e-16\n\njulia> planckreduced(CGS2019) # erg⋅s⋅rad⁻¹\n1.0545718176461565e-27\n\njulia> lightspeed(CGS2019) # cm⋅s⁻¹\n2.99792458e10\n\njulia> permeability(CGS2019) # erg⋅A⁻²⋅cm⁻¹\n12.56637062121048\n\njulia> electronmass(CGS2019 # g\n9.109383701558256e-28\n\n\n\n\n\n","category":"constant"},{"location":"units.html#","page":"Unit systems","title":"Unit systems","text":"Historically, the josephson and klitzing constants have been used to define Conventional and CODATA derived UnitSystem variants.","category":"page"},{"location":"units.html#","page":"Unit systems","title":"Unit systems","text":"julia> josephson(Conventional)\n4.835979e14\n\njulia> klitzing(Conventional)\n25812.807\n\njulia> josephson(CODATA)\n4.835978525e14\n\njulia> klitzing(CODATA)\n25812.8074555","category":"page"},{"location":"units.html#","page":"Unit systems","title":"Unit systems","text":"Geophysics.Conventional\nGeophysics.CODATA","category":"page"},{"location":"units.html#Geophysics.Conventional","page":"Unit systems","title":"Geophysics.Conventional","text":"Conventional\n\nConventional electronic UnitSystem with 1990 tuned josephson and klitzing constants.\n\njulia> boltzmann(Conventional) # J⋅K⁻¹\n1.380649e-23\n\njulia> planckreduced(Conventional) # J⋅s⋅rad⁻¹\n1.054571611438857e-34\n\njulia> lightspeed(Conventional) # m⋅s⁻¹\n2.99792458e8\n\njulia> permeability(Conventional) # H⋅m⁻¹\n1.2566370397608662e-6\n\njulia> electronmass(Conventional) # kg\n9.109383701558256e-31\n\n\n\n\n\n","category":"constant"},{"location":"units.html#Geophysics.CODATA","page":"Unit systems","title":"Geophysics.CODATA","text":"CODATA\n\nMetric UnitSystem based on Committee on Data of the International Science Council.\n\njulia> boltzmann(CODATA) # J⋅K⁻¹\n1.380649e-23\n\njulia> planckreduced(CODATA) # J⋅s⋅rad⁻¹\n1.0545717999940896e-34\n\njulia> lightspeed(CODATA) # m⋅s⁻¹\n2.99792458e8\n\njulia> permeability(CODATA) # H⋅m⁻¹\n1.2566370619358342e-6\n\njulia> electronmass(CODATA) # kg\n9.109383701558256e-31\n\n\n\n\n\n","category":"constant"},{"location":"units.html#","page":"Unit systems","title":"Unit systems","text":"In Britain and the United States an English system of engineering units was commonly used.","category":"page"},{"location":"units.html#","page":"Unit systems","title":"Unit systems","text":"Geophysics.English","category":"page"},{"location":"units.html#Geophysics.English","page":"Unit systems","title":"Geophysics.English","text":"English\n\nEngineering UnitSystem historically used by Britain and United States.\n\njulia> boltzmann(English) # ft⋅lb⋅°R⁻¹\n5.657302466e-24\n\njulia> planckreduced(English) # ft⋅lb⋅s⋅rad⁻¹\n7.778098449204079e-35\n\njulia> lightspeed(English) # ft⋅s⁻¹\n9.83571056e8\n\njulia> permeability(English) # slug⋅ft²⋅?⁻²\n1.2566370614359173e-6\n\njulia> electronmass(English) # kg\n6.241910570285988e-32\n\n\n\n\n\n","category":"constant"},{"location":"units.html#Natural-units-1","page":"Unit systems","title":"Natural units","text":"","category":"section"},{"location":"units.html#","page":"Unit systems","title":"Unit systems","text":"With the introduction of the planckmass a set of natural atomic unit systems can be derived in terms of the gravitational coupling constant.","category":"page"},{"location":"units.html#","page":"Unit systems","title":"Unit systems","text":"julia> αG # (mₑ/mP)^2\n1.751809945750515e-45","category":"page"},{"location":"units.html#","page":"Unit systems","title":"Unit systems","text":"Some of the notable variants include","category":"page"},{"location":"units.html#","page":"Unit systems","title":"Unit systems","text":"Planck       ::UnitSystem{1,1,1,1,√(4π*αG)}\nPlanckGauss  ::UnitSystem{1,1,1,4π,√αG}\nStoney       ::UnitSystem{1,αinv,1,4π,√(αG*αinv)}\nHartree      ::UnitSystem{1,1,αinv,4π/αinv^2,1}\nRydberg      ::UnitSystem{1,1,2αinv,π/αinv^2,1/2}\nSchrodinger  ::UnitSystem{1,1,αinv,4π/αinv^2,√αinv*mₑ/mP}\nElectronic   ::UnitSystem{1,αinv,1,4π,1}\nNatural      ::UnitSystem{1,1,1,1,1}\nNaturalGauss ::UnitSystem{1,1,1,4π,1}\nQCD          ::UnitSystem{1,1,1,1,1/μₚₑ}\nQCDGauss     ::UnitSystem{1,1,1,4π,1/μₚₑ}\nQCDoriginal  ::UnitSystem{1,1,1,4π/αinv,1/μₚₑ}","category":"page"},{"location":"units.html#","page":"Unit systems","title":"Unit systems","text":"Geophysics.Planck\nGeophysics.PlanckGauss\nGeophysics.Stoney\nGeophysics.Hartree\nGeophysics.Rydberg\nGeophysics.Schrodinger\nGeophysics.Electronic\nGeophysics.Natural\nGeophysics.NaturalGauss\nGeophysics.QCD\nGeophysics.QCDGauss\nGeophysics.QCDoriginal","category":"page"},{"location":"units.html#Geophysics.Planck","page":"Unit systems","title":"Geophysics.Planck","text":"Planck\n\nPlanck UnitSystem with the electronmass value √(4π*αG) using gravitational coupling.\n\njulia> boltzmann(Planck)\n1\n\njulia> planckreduced(Planck)\n1\n\njulia> lightspeed(Planck)\n1\n\njulia> permeability(Planck)\n1\n\njulia> electronmass(Planck)\n1.4837079572551133e-22\n\n\n\n\n\n","category":"constant"},{"location":"units.html#Geophysics.PlanckGauss","page":"Unit systems","title":"Geophysics.PlanckGauss","text":"PlanckGauss\n\nPlanck (Gauss) UnitSystem with permeability of 4π and electronmass coupling √αG.\n\njulia> boltzmann(PlanckGauss)\n1\n\njulia> planckreduced(PlanckGauss)\n1\n\njulia> lightspeed(PlanckGauss)\n1\n\njulia> permeability(PlanckGauss)\n12.566370614359172\n\njulia> electronmass(PlanckGauss)\n4.1854628725512725e-23\n\n\n\n\n\n","category":"constant"},{"location":"units.html#Geophysics.Stoney","page":"Unit systems","title":"Geophysics.Stoney","text":"Stoney\n\nStoney UnitSystem with permeability of 4π and electronmass coupling √(αG*αinv).\n\njulia> boltzmann(Stoney)\n1\n\njulia> planckreduced(Stoney)\n137.035999084\n\njulia> lightspeed(Stoney)\n1\n\njulia> permeability(Stoney)\n12.566370614359172\n\njulia> electronmass(Stoney)\n4.899602291219254e-22\n\n\n\n\n\n","category":"constant"},{"location":"units.html#Geophysics.Hartree","page":"Unit systems","title":"Geophysics.Hartree","text":"Hartree\n\nHartree atomic UnitSystem with lightspeed of αinv and permeability of 4π/αinv^2.\n\njulia> boltzmann(Hartree)\n1\n\njulia> planckreduced(Hartree)\n1\n\njulia> lightspeed(Hartree)\n137.035999084\n\njulia> permeability(Hartree)\n0.0006691762566203904\n\njulia> electronmass(Hartree)\n1\n\n\n\n\n\n","category":"constant"},{"location":"units.html#Geophysics.Rydberg","page":"Unit systems","title":"Geophysics.Rydberg","text":"Rydberg\n\nRydberg UnitSystem with lightspeed of 2αinv and permeability of π/αinv^2.\n\njulia> boltzmann(Rydberg)\n1\n\njulia> planckreduced(Rydberg)\n1\n\njulia> lightspeed(Rydberg)\n274.071998168\n\njulia> permeability(Rydberg)\n0.0001672940641550976\n\njulia> electronmass(Rydberg)\n0.5\n\n\n\n\n\n","category":"constant"},{"location":"units.html#Geophysics.Schrodinger","page":"Unit systems","title":"Geophysics.Schrodinger","text":"Schrodinger\n\nSchrodinger UnitSystem with permeability of 4π/αinv^2 and electronmass of √αinv*mₑ/mP.\n\njulia> boltzmann(Schrodinger)\n1\n\njulia> planckreduced(Schrodinger)\n1\n\njulia> lightspeed(Schrodinger)\n137.035999084\n\njulia> permeability(Schrodinger)\n0.0006691762566203904\n\njulia> electronmass(Schrodinger)\n4.899602291219254e-22\n\n\n\n\n\n","category":"constant"},{"location":"units.html#Geophysics.Electronic","page":"Unit systems","title":"Geophysics.Electronic","text":"Electronic\n\nElectronic UnitSystem with planckreduced of αinv and permeability of 4π.\n\njulia> boltzmann(Electronic)\n1\n\njulia> planckreduced(Electronic)\n137.035999084\n\njulia> lightspeed(Electronic)\n1\n\njulia> permeability(Electronic)\n12.566370614359172\n\njulia> electronmass(Electronic)\n1\n\n\n\n\n\n","category":"constant"},{"location":"units.html#Geophysics.Natural","page":"Unit systems","title":"Geophysics.Natural","text":"Natural\n\nNatural UnitSystem with all primary constants having unit value.\n\njulia> boltzmann(Natural)\n1\n\njulia> planckreduced(Natural)\n1\n\njulia> lightspeed(Natural)\n1\n\njulia> permeability(Natural)\n1\n\njulia> electronmass(Natural)\n1\n\n\n\n\n\n","category":"constant"},{"location":"units.html#Geophysics.NaturalGauss","page":"Unit systems","title":"Geophysics.NaturalGauss","text":"NaturalGauss\n\nNatural (Gauss) UnitSystem with the Gaussian permeability value of 4π.\n\njulia> boltzmann(NaturalGauss)\n1\n\njulia> planckreduced(NaturalGauss)\n1\n\njulia> lightspeed(NaturalGauss)\n1\n\njulia> permeability(NaturalGauss)\n12.566370614359172\n\njulia> electronmass(NaturalGauss)\n1\n\n\n\n\n\n","category":"constant"},{"location":"units.html#Geophysics.QCD","page":"Unit systems","title":"Geophysics.QCD","text":"QCD\n\nQunatum chromodynamics UnitSystem with electronmass of 1/μₚₑ or 1/1836.152673432705.\n\njulia> boltzmann(QCD)\n1\n\njulia> planckreduced(QCD)\n1\n\njulia> lightspeed(QCD)\n1\n\njulia> permeability(QCD)\n1\n\njulia> electronmass(QCD)\n0.0005446170214868301\n\n\n\n\n\n","category":"constant"},{"location":"units.html#Geophysics.QCDGauss","page":"Unit systems","title":"Geophysics.QCDGauss","text":"QCDGauss\n\nQunatum chromodynamics (Gauss) UnitSystem with electronmass of 1/μₚₑ.\n\njulia> boltzmann(QCDGauss)\n1\n\njulia> planckreduced(QCDGauss)\n1\n\njulia> lightspeed(QCDGauss)\n1\n\njulia> permeability(QCDGauss)\n12.566370614359172\n\njulia> electronmass(QCDGauss)\n0.0005446170214868301\n\n\n\n\n\n","category":"constant"},{"location":"units.html#Geophysics.QCDoriginal","page":"Unit systems","title":"Geophysics.QCDoriginal","text":"QCDoriginal\n\nQunatum chromodynamics (original) UnitSystem with permeability of 4π/αinv.\n\njulia> boltzmann(QCDoriginal)\n1\n\njulia> planckreduced(QCDoriginal)\n1\n\njulia> lightspeed(QCDoriginal)\n1\n\njulia> permeability(QCDoriginal)\n0.09170123688926636\n\njulia> electronmass(QCDoriginal)\n0.0005446170214868301\n\n\n\n\n\n","category":"constant"},{"location":"units.html#Fundamental-constants-of-physics-1","page":"Unit systems","title":"Fundamental constants of physics","text":"","category":"section"},{"location":"units.html#","page":"Unit systems","title":"Unit systems","text":"The following are fundamental constants of physics:","category":"page"},{"location":"units.html#","page":"Unit systems","title":"Unit systems","text":"avogadro\nboltzmann\nuniversal\nplanck\nplanckreduced\nlightspeed\nplanckmass\natomicmass\nprotonmass\nelectronmass\nnewton\neinstein\npermeability\npermittivity\ncoulomb\nstefan\nradiationintensity\nimpedance\ncharge\nmagneton\nconductance\nfaraday\nmagneticflux\njosephson\nklitzing\nhardtree\nrydberg\nbohr\nbohrreduced","category":"page"},{"location":"units.html#Geophysics.avogadro","page":"Unit systems","title":"Geophysics.avogadro","text":"avogadro(x) = universal(x)/boltzmann(x) # Mᵤ/atomicmass(x), Mᵤ ≈ 0.001-3.5e-13\n\nAvogadro NA is molarmass(x)/molecularmass(x) number of atoms in a 12 g sample of C₁₂.\n\njulia> avogadro(Metric) # mol⁻¹\n6.02214076e23\n\njulia> avogadro(English) # slug-mol⁻¹\n8.788653773538476e24\n\n\n\n\n\n","category":"function"},{"location":"units.html#Geophysics.boltzmann","page":"Unit systems","title":"Geophysics.boltzmann","text":"boltzmann(x) = universal(x)/avogadro(x)\n\nBoltzmann constant kB is the entropy amount of a unit number microstate permutation.\n\npressure*molecularmass == density*boltzmann*temperature\n\nIt satisfies the ideal gas law.\n\njulia> boltzmann(Metric) # J⋅K⁻¹\n1.380649e-23\n\njulia> boltzmann(English) # ft⋅lb⋅°R⁻¹\n5.657302466e-24\n\n\n\n\n\n","category":"function"},{"location":"units.html#Geophysics.universal","page":"Unit systems","title":"Geophysics.universal","text":"universal(x) = boltzmann(x)*avogadro(x)\n\nUniversal gas constant Rᵤ is factored into specific gasconstant(x)*molarmass(x) values.\n\npressure*molarmass == density*universal*temperature\n\nIt satisfies the ideal gas law.\n\njulia> universal(Metric) # J⋅K⁻¹⋅mol⁻¹\n8.31446261815324\n\njulia> universal(English) # ft⋅lb⋅°R⁻¹⋅slug-mol⁻¹\n49.720072665859426\n\n\n\n\n\n","category":"function"},{"location":"units.html#Geophysics.planck","page":"Unit systems","title":"Geophysics.planck","text":"planck(x) = 2π*planckreduced(x)\n\nPlanck constant 𝘩 is energy per electromagnetic frequency (J⋅s or ft⋅lb⋅s).\n\njulia> planck(Metric) # J⋅s\n6.62607015e-34\n\njulia> planck(English) # ft⋅lb⋅s\n4.88712338938354e-34\n\n\n\n\n\n","category":"function"},{"location":"units.html#Geophysics.planckreduced","page":"Unit systems","title":"Geophysics.planckreduced","text":"planckreduced(x) = planck(x)/2π\n\nReduced Planck constant ħ is a Planck per radian (J⋅s⋅rad⁻¹ or ft⋅lb⋅s⋅rad⁻¹).\n\njulia> planckreduced(Metric) # J⋅s⋅rad⁻¹\n1.0545718176461565e-34\n\njulia> planckreduced(English) # ft⋅lb⋅s⋅rad⁻¹\n7.778098449204079e-35\n\n\n\n\n\n","category":"function"},{"location":"units.html#Geophysics.lightspeed","page":"Unit systems","title":"Geophysics.lightspeed","text":"lightspeed(x) = 1/sqrt(μ₀*ε₀)\n\nSpeed of light in a vacuum 𝘤 for massless particles (m⋅s⁻¹ or ft⋅s⁻¹).\n\njulia> universal(Metric) # m⋅s⁻¹\n2.99792458e8\n\njulia> universal(English) # ft⋅s⁻¹\n9.83571056e8\n\n\n\n\n\n","category":"function"},{"location":"units.html#Geophysics.planckmass","page":"Unit systems","title":"Geophysics.planckmass","text":"planckmass(U::UnitSystem) = mass(2.176434e-8,U)\n\nPlanck mass factor mP from the gravitational coupling constant αG (kg or slugs).\n\njuila> newton(Metric) # m³⋅kg⁻¹⋅s⁻²\n6.674302101972536e-11\n\njulia> newton(English) # ft³⋅slug⁻¹⋅s⁻²\n3.439783265696398e-8\n\n\n\n\n\n","category":"function"},{"location":"units.html#Geophysics.atomicmass","page":"Unit systems","title":"Geophysics.atomicmass","text":"atomicmass(U::UnitSystem) = Mᵤ/avogadro(U) # 0.000999999999656256 ≈ 0.001-3.5e-13\n\nAtomic mass unit mᵤ of 1/12 of the C₁₂ carbon-12 atom's mass  (kg or slugs).\n\njulia> atomicmass(Metric) # kg\n1.6605390666030467e-27\n\njulia> atomicmass(English) # slugs\n1.137830691052058e-28\n\n\n\n\n\n","category":"function"},{"location":"units.html#Geophysics.protonmass","page":"Unit systems","title":"Geophysics.protonmass","text":"protonmass(U::UnitSystem) = 1.007276466621atomicmass(U)\n\nProton mass unit mₚ of subatomic particle with +𝘦 elementary charge  (kg or slugs).\n\njulia> protonmass(Metric) # kg\n1.6726219236940502e-27\n\njulia> protonmass(English) # slugs\n1.1461100780958476e-28\n\n\n\n\n\n","category":"function"},{"location":"units.html#Geophysics.electronmass","page":"Unit systems","title":"Geophysics.electronmass","text":"electronmass(U::UnitSystem) = protonmass(U)/1836.152673432705\n\nElectron rest mass unit mₑ of subatomic particle with -𝘦 elementary charge  (kg or slugs).\n\njulia> electronmass(Metric) # kg\n9.109383701558256e-31\n\njulia> electronmass(English) # slugs\n6.241910570285988e-32\n\n\n\n\n\n","category":"function"},{"location":"units.html#Geophysics.newton","page":"Unit systems","title":"Geophysics.newton","text":"newton(x) = lightspeed(x)*planckreduced(x)/planckmass(x)^2\n\nUniversal gravitational constant GG of Newton's law (m³⋅kg⁻¹⋅s⁻² or ft³⋅slug⁻¹⋅s⁻²).\n\njuila> newton(Metric) # m³⋅kg⁻¹⋅s⁻²\n6.674302101972536e-11\n\njulia> newton(English) # ft³⋅slug⁻¹⋅s⁻²\n3.439783265696398e-8\n\n\n\n\n\n","category":"function"},{"location":"units.html#Geophysics.einstein","page":"Unit systems","title":"Geophysics.einstein","text":"einstein(U::UnitSystem) = 8π*newton(U)/lightspeed(U)^4\n\nEinstein's gravitational constant from the Einstein field equations (? or ?).\n\njulia> einstein(Metric) # ?\n1.6605390666030467e-27\n\njulia> einstein(English) # ?\n1.137830691052058e-28\n\n\n\n\n\n","category":"function"},{"location":"units.html#Geophysics.permeability","page":"Unit systems","title":"Geophysics.permeability","text":"permeability(x) = 4π*1e-7\n\nMagnetic permeability in a classical vacuum defined as μ₀ in SI units (H⋅m⁻¹).\n\njulia> permeability(Metric) # H⋅m⁻¹\n1.2566370614359173e-6\n\njulia> permeability(English) # slug⋅ft²⋅?⁻²\n1.2566370614359173e-6\n\n\n\n\n\n","category":"function"},{"location":"units.html#Geophysics.permittivity","page":"Unit systems","title":"Geophysics.permittivity","text":"permittivity(x) = 1/μ₀/lightspeed(x)^2\n\nDielectric permittivity constant ε₀ of a classical vacuum (C²⋅N⁻¹⋅m⁻²).\n\njulia> permittivity(Metric) # C²⋅N⁻¹⋅m⁻²\n8.854187817620389e-12\n\n\n\n\n\n","category":"function"},{"location":"units.html#Geophysics.coulomb","page":"Unit systems","title":"Geophysics.coulomb","text":"coulomb(x) = 1/4π/ϵ₀\n\nElectrostatic proportionality constant kₑ for the Coulomb's law force (N⋅m²⋅C⁻²).\n\njulia> coulomb(Metric) # N⋅m²⋅C⁻²\n8.987551787368177e9\n\n\n\n\n\n","category":"function"},{"location":"units.html#Geophysics.stefan","page":"Unit systems","title":"Geophysics.stefan","text":"stefan(x) = 2π^5*boltzmann(U)^4/(15planck(U)^3*lightspeed(U)^2)\n\nStefan-Boltzmann proportionality σ of black body radiation (W⋅m⁻²⋅K⁻⁴ or ?⋅ft⁻²⋅°R⁻⁴).\n\njulia> stefan(Metric) # W⋅m⁻²⋅K⁻⁴\n5.6703744191844314e-8\n\njulia> stefan(English) # lb⋅s⁻¹⋅ft⁻³⋅°R⁻⁴\n3.701300130852963e-10\n\n\n\n\n\n","category":"function"},{"location":"units.html#Geophysics.impedance","page":"Unit systems","title":"Geophysics.impedance","text":"impedance(U::UnitSystem) = permeability(U)*lightspeed(U)\n\nVacuum impedance of free space Z₀ is magnitude ratio of electric to magnetic field (Ω).\n\njulia> impedance(Metric) # Ω\n376.73031346177066\n\n\n\n\n\n","category":"function"},{"location":"units.html#Geophysics.charge","page":"Unit systems","title":"Geophysics.charge","text":"charge(U::UnitSystem) = sqrt(2𝘩/137.035999084impedance(U))\n\nQuantized elementary charge 𝘦 of a proton or electron  (C).\n\njulia> charge(Metric) # C\n1.6021766344367608e-19\n\n\n\n\n\n","category":"function"},{"location":"units.html#Geophysics.magneton","page":"Unit systems","title":"Geophysics.magneton","text":"magneton(U::UnitSystem) = charge(U)*planckreduced(U)/2electronmass(U)\n\nBohr magneton μB natural unit for expressing magnetic moment of electron (J⋅T⁻¹).\n\njulia> magneton(Metric) # J⋅T⁻¹\n9.274010080830994e-24\n\n\n\n\n\n","category":"function"},{"location":"units.html#Geophysics.conductance","page":"Unit systems","title":"Geophysics.conductance","text":"conductance(U::UnitSystem) = 2charge(U)^2/𝘩\n\nConductance quantum G₀ is a quantized unit of electrical conductance (S).\n\njulia> conductance(Metric) # S\n7.748091734087984e-5\n\n\n\n\n\n","category":"function"},{"location":"units.html#Geophysics.faraday","page":"Unit systems","title":"Geophysics.faraday","text":"faraday(U::UnitSystem) = charge(U)*avogadro(U)\n\nElectric charge per mole of electrons 𝔉 based on elementary charge (C⋅mol⁻¹).\n\njulia> faraday(Metric) # C⋅mol⁻¹\n96485.33214961237\n\n\n\n\n\n","category":"function"},{"location":"units.html#Geophysics.magneticflux","page":"Unit systems","title":"Geophysics.magneticflux","text":"magneticflux(U::UnitSystem) = planck(U)/2charge(U)\n\nMagnetic flux quantum Φ₀ is 1/josephson(U) (Wb).\n\njulia> magneticflux(Metric) # Wb\n2.067833847898228e-15\n\n\n\n\n\n","category":"function"},{"location":"units.html#Geophysics.josephson","page":"Unit systems","title":"Geophysics.josephson","text":"josephson(U::UnitSystem) = 2charge(U)/planck(U)\n\nJosephson constant KJ relating potential difference to irradiation frequency across junction (Hz⋅V⁻¹).\n\njulia> josephson(Metric) # Hz⋅V⁻¹\n4.835978485488147e14\n\n\n\n\n\n","category":"function"},{"location":"units.html#Geophysics.klitzing","page":"Unit systems","title":"Geophysics.klitzing","text":"klitzing(U::UnitSystem) = 2/conductance(U)\n\nQuantized Hall resistance RK (Ω).\n\njulia> klitzing(Metric) # Ω\n25812.80744523112\n\n\n\n\n\n","category":"function"},{"location":"units.html#Geophysics.hardtree","page":"Unit systems","title":"Geophysics.hardtree","text":"hardtree(U::UnitSystem) = electronmass(U)*(lightspeed(U)/137.035999084)^2\n\nHardtree electric potential energy Eₕ of the hydrogen atom at ground state (J).\n\njulia> hardtree(Metric) # J\n4.359744722207209e-18\n\n\n\n\n\n","category":"function"},{"location":"units.html#Geophysics.rydberg","page":"Unit systems","title":"Geophysics.rydberg","text":"rydberg(U::UnitSystem) = hardtree(U)/2planck(U)/lightspeed(U)\n\nRydberg constant R∞ is lowest energy photon capable of ionizing atom at ground state (m⁻¹).\n\njulia> rydberg(Metric) # m⁻¹\n1.09737315681601e7\n\n\n\n\n\n","category":"function"},{"location":"units.html#Geophysics.bohr","page":"Unit systems","title":"Geophysics.bohr","text":"bohr(U) = 137.035999084*planckreduced(U)/electronmass(U)/lightspeed(U)\n\nBohr radius of the hydrogen atom in its ground state a₀ (m).\n\njulia> bohr(Metric) # m\n5.291772109022829e-11\n\n\n\n\n\n","category":"function"},{"location":"units.html#Geophysics.bohrreduced","page":"Unit systems","title":"Geophysics.bohrreduced","text":"bohrreduced(U::UnitSystem) = electronmass(U)/bohr(U)/1836.152673432705\n\nReduced Bohr radius including the effect of reduced mass in hydrogen atom (m).\n\njulia> bohrreduced(Metric) # m\n2.6253145122281045e-44\n\n\n\n\n\n","category":"function"},{"location":"units.html#Common-conversion-factors-1","page":"Unit systems","title":"Common conversion factors","text":"","category":"section"},{"location":"units.html#","page":"Unit systems","title":"Unit systems","text":"Common conversion factors include moles, molecules, kilograms, slugs, meters, feet, kelvin, and rankine.","category":"page"},{"location":"units.html#","page":"Unit systems","title":"Unit systems","text":"moles\nmolecules\nkilograms\nslugs\nmeters\nfeet\nkelvin\nrankine","category":"page"},{"location":"units.html#Geophysics.moles","page":"Unit systems","title":"Geophysics.moles","text":"moles(N::Real,U::UnitSystem=Metric) = N/avogadro(U)\n\nConverts the number of molecules N to number of moles (mol).\n\n\n\n\n\n","category":"function"},{"location":"units.html#Geophysics.molecules","page":"Unit systems","title":"Geophysics.molecules","text":"molecules(n::Real,U::UnitSystem=Metric) = n*avogadro(U)\n\nConverts the number of moles n to number of molecules (dimensionless).\n\n\n\n\n\n","category":"function"},{"location":"units.html#Geophysics.kilograms","page":"Unit systems","title":"Geophysics.kilograms","text":"kilograms(m::Real) = 14.593902938825488m\n\nConverts mass m from slugs to kilogram (kg).\n\n\n\n\n\n","category":"function"},{"location":"units.html#Geophysics.slugs","page":"Unit systems","title":"Geophysics.slugs","text":"slugs(m::Real) = 0.06852176584918959m\n\nConverts mass m from kilograms to slugs (slug).\n\n\n\n\n\n","category":"function"},{"location":"units.html#Geophysics.meters","page":"Unit systems","title":"Geophysics.meters","text":"meters(d) = 0.3048000001333915d\n\nConverts distance d from feet to meters (m).\n\n\n\n\n\n","category":"function"},{"location":"units.html#Geophysics.feet","page":"Unit systems","title":"Geophysics.feet","text":"feet(d) = 3.2808398935773093d\n\nConverts distance d from meters to feet (ft).\n\n\n\n\n\n","category":"function"},{"location":"units.html#Geophysics.kelvin","page":"Unit systems","title":"Geophysics.kelvin","text":"kelvin(T) = (5/9)T\n\nConverts temperature T from degrees Rankine to Kelvin (K).\n\n\n\n\n\n","category":"function"},{"location":"units.html#Geophysics.rankine","page":"Unit systems","title":"Geophysics.rankine","text":"rankine(T) = (9/5)T\n\nConverts temperature T from Kelvin to degrees Rankine (°R).\n\n\n\n\n\n","category":"function"},{"location":"units.html#Index-1","page":"Unit systems","title":"Index","text":"","category":"section"},{"location":"units.html#","page":"Unit systems","title":"Unit systems","text":"Pages = [\"units.md\"]","category":"page"},{"location":"index.html#Geophysics.jl-1","page":"Home","title":"Geophysics.jl","text":"","category":"section"},{"location":"index.html#","page":"Home","title":"Home","text":"Planetary science data for atmospheric geophysical models","category":"page"},{"location":"index.html#","page":"Home","title":"Home","text":"(Image: DOI) (Image: Build Status) (Image: Build status) (Image: Coverage Status) (Image: codecov.io) (Image: Gitter)","category":"page"},{"location":"index.html#","page":"Home","title":"Home","text":"Provides Atmosphere models based on Air Research and Development Command ARDC and the United States (1922, 1925, 1956, 1959, 1962, 1966, 1976) Standard Atmosphere US22,US25,US56,US59,US62,US66,US76 available also in English units US22E,US25E,US56E,US59E,US62E,US66E,US76E. Provided the local absolute sea level and gravitational acceleration, the Weather can be initialized based on temperature and pressure. Presets for the Standard atmosphere are provided: Earth1922, Earth1925, Earth1956, Earth1959, Earth1962, Earth1966, Earth1976, Earth1922English, Earth1925English, Earth1956English, Earth1959English, Earth1962English, Earth1966English, Earth1976English. By default the 1959 model with metric units is used for Standard atmosphere, although a different year can be specified with environment variable STDATM and the default unit system can be specified with the GEOUNITS environment variable.","category":"page"},{"location":"index.html#","page":"Home","title":"Home","text":"julia> using Geophysics\n\njulia> h = 1000 # altitude, m\n1000\n\njulia> gravity(h)\n9.803565306802405\n\njulia> temperature(h)\n281.66102237169474\n\njulia> pressure(h)\n89876.28158431675\n\njulia> sonicspeed(h)\n336.4347118683662","category":"page"},{"location":"index.html#","page":"Home","title":"Home","text":"Values which can be obtained at geometric altitude include gravity, temperature, pressure, density, sonicspeed, conductivity, viscosity, kinematic, volume, energy, enthalpy, heatcapacity, diffusivity, prandtl, and impedance. In the future, more varieties of atmospheric models will be added for various planets along with winds aloft and turbulent gust distribution data. Weather data from internet sources may be imported in the future.","category":"page"},{"location":"index.html#","page":"Home","title":"Home","text":"Pages = [\"units.md\",\"references.md\"]","category":"page"},{"location":"index.html#","page":"Home","title":"Home","text":"This package is not limited to atmospheric data: other geophysical data features are intended to be added for oceans, temperature and pressure inside the planets, as well as electrical and magnetic properties of planets. In this package, any simple Geophysical properties of planets may be added. Other simple geophysical data about planets, can be added in a collaborative effort. Complicated models will be excluded from this package, as it is only intended to provide a minimal foundation for geophysical data and constants of various planets, more complicated models should be built separately in packages to build on Geophysics. For example, some geographic conditions can be calculated externally, and then Geophysics is used to load that data.","category":"page"},{"location":"index.html#","page":"Home","title":"Home","text":"This Geophysics package for the Julia language was created by github.com/chakravala for research.","category":"page"}]
}
