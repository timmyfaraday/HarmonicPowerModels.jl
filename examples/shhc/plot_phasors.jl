################################################################################
# Companion plot: deterministic-equivalence phasor diagrams (Argand plane).       #
# Overlays the dHHC / SOC / QCQP current and voltage phasors per harmonic — if    #
# the three formulations are equivalent the phasor tips coincide.                 #
#                                                                               #
# Requires Plots.jl (not a HarmonicPowerModels dependency):                       #
#   julia> ENV["SHHC_NET"]="industrial"; include("plot_phasors.jl")             #
# Produces output/<net>_phasors_V.pdf and _phasors_I.pdf.                        #
################################################################################
using Plots, DelimitedFiles, Printf

const NET = get(ENV, "SHHC_NET", "industrial")
outdir = joinpath(@__DIR__, "output")
raw, _ = readdlm(joinpath(outdir, "$(NET)_det_phasors.csv"), ','; header = true)
# columns: method, quantity, id, h, re, im, mag, ang_deg
methods = ["dHHC", "SOC", "QCQP"]
mstyle  = Dict("dHHC"=>(:circle,:steelblue,8), "SOC"=>(:diamond,:orange,6), "QCQP"=>(:xcross,:black,6))

function plot_quantity(q, fname, title)
    harms = sort(unique(Int.(raw[raw[:,2].==q, 4])))
    plts = []
    for h in harms
        p = plot(title = "h = $h", legend = (h == harms[1]), aspect_ratio = :equal,
                 xlabel = "real", ylabel = "imag", framestyle = :origin)
        for meth in methods
            idx = (raw[:,1] .== meth) .& (raw[:,2] .== q) .& (Int.(raw[:,4]) .== h)
            re = Float64.(raw[idx, 5]); im = Float64.(raw[idx, 6])
            for k in eachindex(re)                      # phasor arrows from origin
                plot!(p, [0, re[k]], [0, im[k]]; color = mstyle[meth][2], lw = 1, label = "")
            end
            mk, col, ms = mstyle[meth]
            scatter!(p, re, im; marker = mk, color = col, ms = ms, label = meth, markerstrokewidth = 0)
        end
        push!(plts, p)
    end
    fig = plot(plts...; layout = (2, 2), size = (1000, 900), plot_title = title)
    savefig(fig, joinpath(outdir, fname))
    println("✓ wrote ", fname)
end

gr()
plot_quantity("V", "$(NET)_phasors_V.pdf", "Bus voltage phasors @E[φ] — $(uppercase(NET))")
plot_quantity("I", "$(NET)_phasors_I.pdf", "Injected current phasors @E[φ] — $(uppercase(NET))")
