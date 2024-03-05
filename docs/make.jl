using Documenter
using HarmonicPowerModels

makedocs(
    modules     = [HarmonicPowerModels],
    format      = Documenter.HTML(mathengine = Documenter.MathJax()),
    sitename    = "HarmonicPowerModels.jl",
    authors     = "Tom Van Acker",
    pages       =   [    "Home"                   => "index.md",
                         "Components"             => 
                         [    "Bus"               => "bus.md",
                              "Reference Bus"     => "ref_bus.md",
                              "Branch"            => "branch.md",
                              "Transformer"       => "xfmr.md",
                              "Harmonic Load"     => "load.md"
                         ],
                         "Problem Formulation"    =>
                         [    "HPF"               => "power_flow.md",
                              "HOPF"              => "optimal_power_flow.md",
                              "HHC"               => "hosting_capacity.md"
                         ]
                    ]
)

deploydocs(
     repo = "github.com/timmyfaraday/HarmonicPowerModels.jl.git"
)