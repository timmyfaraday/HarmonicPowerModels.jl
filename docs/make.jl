push!(LOAD_PATH,"../src/")
using Documenter, HarmonicPowerModels

# A flag to check if we are running in a GitHub action.
const _IS_GITHUB_ACTIONS = get(ENV, "GITHUB_ACTIONS", "false") == "true"

# Pass --pdf to build the PDF. On GitHub actions, we always build the PDF.
const _PDF = findfirst(isequal("--pdf"), ARGS) !== nothing || _IS_GITHUB_ACTIONS

const _PAGES = [
    "Home"          => "index.md",
    "Transformer"   => "xfmr.md"
];

@time Documenter.makedocs(
    modules = [HarmonicPowerModels],
    format = Documenter.HTML(mathengine = Documenter.MathJax()),
    sitename = "HarmonicPowerModels",
    authors = "Tom Van Acker, Frederik Geth and contributors.",
    pages = _PAGES,
    )

if _PDF
    # latex_platform = _IS_GITHUB_ACTIONS ? "docker" : "native"
    latex_platform = "docker"
    @time Documenter.makedocs(
        sitename = "HarmonicPowerModels",
        authors = "The HarmonicPowerModels core developers and contributors",
        format = Documenter.LaTeX(platform = latex_platform),
        build = "latex_build",
    )
    # Hack for deploying: copy the pdf (and only the PDF) into the HTML build
    # directory! We don't want to copy everything in `latex_build` because it
    # includes lots of extraneous LaTeX files.
    cp(
        joinpath(@__DIR__, "latex_build", "HarmonicPowerModels.pdf"),
        joinpath(@__DIR__, "build", "HarmonicPowerModels.pdf"),
    )
end

Documenter.deploydocs(
    repo = "github.com/timmyfaraday/HarmonicPowerModels.jl.git",
    push_preview = true,
)