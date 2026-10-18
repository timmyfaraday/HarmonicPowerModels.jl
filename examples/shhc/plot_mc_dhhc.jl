################################################################################
# Companion plot for mc_dhhc_validation.jl — renders the MC hosting-capacity      #
# distribution with the reference levels (q05/q95/median, dHHC worst case, sHHC   #
# Cantelli/Gaussian) overlaid as vertical lines.                                  #
#                                                                               #
# Requires Plots.jl in the active environment (NOT a HarmonicPowerModels dep):    #
#     julia> import Pkg; Pkg.add("Plots")                                         #
#     julia> ENV["SHHC_NET"]="industrial"; include("plot_mc_dhhc.jl")            #
# Produces output/mc_dhhc_<net>.pdf (and .png).                                   #
################################################################################
using Plots, DelimitedFiles, Printf

const NET = get(ENV, "SHHC_NET", "industrial")
outdir = joinpath(@__DIR__, "output")

samples, _ = readdlm(joinpath(outdir, "mc_dhhc_$(NET).csv"), ','; header=true)
obj = Float64.(samples[:, 2])

refs, _ = readdlm(joinpath(outdir, "mc_dhhc_$(NET)_refs.csv"), ','; header=true)
R = Dict(string(refs[i,1]) => Float64(refs[i,2]) for i in 1:size(refs,1))

gr()
h = histogram(obj; bins=40, normalize=:probability, legend=:topright,
              fillalpha=0.55, linecolor=:white, color=:steelblue,
              xlabel="Harmonic hosting capacity Σ I_h [pu]",
              ylabel="probability", title="MC dHHC hosting distribution — $(uppercase(NET))",
              size=(900,520))

style = Dict("q05"=>(:black,:solid,3), "median"=>(:gray,:dash,2), "q95"=>(:black,:dot,2),
             "dHHC_Ephi"=>(:firebrick,:solid,2),
             "sHHC_Cantelli"=>(:orange,:dashdot,3), "sHHC_Gaussian"=>(:goldenrod,:dashdot,3),
             "sHHC_empirical_SOC"=>(:seagreen,:solid,3), "sHHC_empirical_QCQP"=>(:purple,:dot,2))
labels = Dict("q05"=>"q05 (95%-reliable)", "median"=>"median", "q95"=>"q95",
              "dHHC_Ephi"=>"dHHC @E[φ]", "sHHC_Cantelli"=>"sHHC Cantelli",
              "sHHC_Gaussian"=>"sHHC Gaussian", "sHHC_empirical_SOC"=>"sHHC empirical λ* (SOC)",
              "sHHC_empirical_QCQP"=>"sHHC empirical λ* (QCQP)")
for k in ["dHHC_Ephi","q05","median","q95","sHHC_Cantelli","sHHC_Gaussian",
          "sHHC_empirical_SOC","sHHC_empirical_QCQP"]
    haskey(R,k) || continue
    c,ls,lw = style[k]
    vline!([R[k]]; color=c, linestyle=ls, linewidth=lw, label=labels[k])
end

savefig(h, joinpath(outdir, "mc_dhhc_$(NET).pdf"))
savefig(h, joinpath(outdir, "mc_dhhc_$(NET).png"))
println("✓ wrote output/mc_dhhc_$(NET).pdf / .png")
