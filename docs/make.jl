using Documenter
using HarmonicPowerModels

makedocs(
    modules     = [HarmonicPowerModels],
    format      = Documenter.HTML(mathengine = Documenter.MathJax()),
    sitename    = "HarmonicPowerModels.jl",
    authors     = "Tom Van Acker",
    pages       =   [ "Home"              => "index.md"
                    ]
)

deploydocs(
     repo = "github.com/timmyfaraday/HarmonicPowerModels.jl.git"
)