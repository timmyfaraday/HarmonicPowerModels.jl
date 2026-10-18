################################################################################
# Companion plot: stochastic-equivalence MC distribution with sHHC overlays.      #
# Histogram of the Monte-Carlo dHHC objective with vertical lines for q05/median/ #
# q95 and the sHHC_SOC / sHHC_QCQP outcomes under each chance bound.              #
#                                                                               #
# Requires Plots.jl:                                                              #
#   julia> ENV["SHHC_NET"]="industrial"; include("plot_mc_overlay.jl")          #
# Produces output/<net>_mc_overlay.pdf / .png.                                   #
################################################################################
using Plots, DelimitedFiles, Printf

const NET = get(ENV, "SHHC_NET", "industrial")
outdir = joinpath(@__DIR__, "output")
mc, _   = readdlm(joinpath(outdir, "$(NET)_stoch_mc.csv"), ','; header = true)
obj     = Float64.(mc[:, 2])
refs, _ = readdlm(joinpath(outdir, "$(NET)_stoch_refs.csv"), ','; header = true)
R = Dict(string(refs[i,1]) => Float64(refs[i,2]) for i in 1:size(refs,1))

gr()
h = histogram(obj; bins = 40, normalize = :probability, color = :steelblue, fillalpha = 0.5,
              linecolor = :white, legend = :topright, size = (1000, 560), label = "MC dHHC",
              xlabel = "Σ_h I_h per unit [pu]", ylabel = "probability",
              title = "MC dHHC objective vs sHHC — $(uppercase(NET))")

# quantile guides
for (k, c, ls) in (("q05",:black,:solid), ("median",:gray,:dash), ("q95",:black,:dot))
    haskey(R,k) && vline!([R[k]]; color=c, linestyle=ls, lw=2, label=k)
end
# deterministic dHHC @E[φ] (no-uncertainty reference)
haskey(R,"dHHC_Ephi") && vline!([R["dHHC_Ephi"]]; color=:purple, linestyle=:solid, lw=2.5, label="dHHC @E[φ]")
# sHHC overlays: SOC solid, QCQP dotted; colour per bound
bcol = Dict("empirical"=>:seagreen, "gaussian"=>:goldenrod, "cantelli"=>:firebrick)
for b in ("empirical","gaussian","cantelli"), (mdl, ls) in (("SOC",:solid), ("QCQP",:dashdot))
    k = "$(mdl)_$b"
    haskey(R,k) && vline!([R[k]]; color=bcol[b], linestyle=ls, lw=2, label="$mdl $b")
end

savefig(h, joinpath(outdir, "$(NET)_mc_overlay.pdf"))
savefig(h, joinpath(outdir, "$(NET)_mc_overlay.png"))
println("✓ wrote $(NET)_mc_overlay.pdf / .png")
