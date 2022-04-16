#   This file is part of Geophysics.jl. It is licensed under the AGPL license
#   Geophysics Copyright (C) 2019 Michael Reed

using Documenter, AbstractTensors, Geophysics, MeasureSystems

makedocs(
    modules = [Geophysics,MeasureSystems,UnitSystems],
    doctest = false,
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    sitename = "Geophysics.jl",
    authors = "Michael Reed",
    pages = Any[
        "Home" => "index.md",
        "UnitSystems.jl" => Any[
            "unitsystems.md",
            "constants.md",
            "convert.md",
            "units.md"
           ],
        "References" => "references.md"
        ]
)

struct MyDeployConfig <: Documenter.DeployConfig end

Documenter.deploy_folder(cfg::MyDeployConfig; repo, devbranch, push_preview, devurl, kwargs...) = devurl

Documenter.authentication_method(::MyDeployConfig) = Documenter.SSH

deploydocs(
    repo   = "github.com/chakravala/Geophysics.jl.git",
    deploy_config = MyDeployConfig(),
    target = "build",
    deps = nothing,
    make = nothing
)
