using Documenter, RouteSequencing

DocMeta.setdocmeta!(RouteSequencing, :DocTestSetup, :(using RouteSequencing); recursive=true)

format_HTML = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true",
                         collapselevel = 1,
                         assets = ["assets/juliareach-light.css"])

format_LATEX = Documenter.LaTeX()

makedocs(
    sitename = "RouteSequencing.jl",
    modules = [RouteSequencing],
    format = format_HTML, # select output format 
    pages = [
        "Introduction" => "index.md",
        "Theory" => Any["Hamiltonian paths" => "manual/hamiltonian.md"],
        "API Reference" => Any["Routes and Stops"=>"lib/routes.md",
                               "Solvers"=>"lib/solvers.md",
                               "Utilities"=>"lib/utils.md"],
        "About" => "about.md"
    ],
    strict = false
)

deploydocs(
    repo = "github.com/mforets/RouteSequencing.jl.git",
    push_preview=true
)
