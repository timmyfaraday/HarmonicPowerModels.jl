################################################################################
# HarmonicPowerModels.jl — SHHC paper §IV-C3: choice of the stochastic bound     #
#                                                                               #
# sHHC_SOC with the Cantelli, Gaussian and empirical risk factor against the     #
# exact (1−ε) quantile, for one and four converter types on the purely          #
# inductive base network, N = 21, and ε = 0.01 … 0.20. For each solution, the  #
# violation probability of the binding individual harmonic voltage limit is    #
# evaluated by Monte Carlo (20 000 samples).                                     #
#                                                                               #
# Output: output/paper/bounds_results.csv                                        #
# Author: Tom Van Acker                                                          #
################################################################################

include(joinpath(@__DIR__, "paper_common.jl"))

data  = ieee30_scenario("inductive", ipopt)
units = paper_units(data)
io = open(joinpath(PAPER_OUT, "bounds_results.csv"), "w")
csvline(io, "T", "eps", "bound", "hosting", "exact", "viol")
for T in (1, 4), eps in (0.01, 0.02, 0.05, 0.10, 0.20)
    types = types_spread(units, T)
    exact = tot(recombination("inductive", NS, "local"; data = data, types = types, prob = 1 - eps), :q95)
    for b in ("cantelli", "gaussian", "empirical")
        I = try soc_hosting(data, units, types, b; eps = eps) catch; nothing end
        I === nothing && (csvline(io, T, eps, b, "NaN", @sprintf("%.6e", exact), "NaN"); continue)
        r = recombination("inductive", NS, "local"; data = data, types = types, eval_I = I)
        csvline(io, T, eps, b, @sprintf("%.6e", sum(values(I))), @sprintf("%.6e", exact),
                @sprintf("%.4f", maximum(r[h][:viol] for h in harm)))
    end
end
close(io)
println("✓ output/paper/bounds_results.csv")
