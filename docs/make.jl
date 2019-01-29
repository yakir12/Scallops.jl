using Documenter, Scallops

makedocs(
    modules = [Scallops],
    format = :html,
    sitename = "Scallops.jl",
    pages = Any["index.md"]
)

deploydocs(
    repo = "github.com/yakir12/Scallops.jl.git",
    target = "build",
    julia = "1.0",
    deps = nothing,
    make = nothing,
)
