################################################################################
# HarmonicPowerModels.jl — SHHC paper §IV-C3: choice of the stochastic bound     #
#                                                                               #
# sHHC_SOC with the Cantelli, Gaussian and empirical risk factor against the     #
# exact (1−ε) quantile, for the diode, thyristor and mixed population on the    #
# purely inductive base network, N = 21, and ε = 0.01 … 0.20. For each solution, the  #
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
csvline(io, "population", "eps", "bound", "hosting", "exact", "viol")
for pop in POPULATIONS, eps in (0.01, 0.02, 0.05, 0.10, 0.20)
    ranges, shapes = population(units, pop)
    exact = tot(recombination("inductive", NS, "local"; data = data, ranges = ranges, shapes = shapes, prob = 1 - eps), :q95)
    for b in ("cantelli", "gaussian", "empirical")
        I = try soc_hosting(data, units, nothing, b; eps = eps, pop = pop) catch; nothing end
        I === nothing && (csvline(io, pop, eps, b, "NaN", @sprintf("%.6e", exact), "NaN"); continue)
        r = recombination("inductive", NS, "local"; data = data, ranges = ranges, shapes = shapes, eval_I = I)
        csvline(io, pop, eps, b, @sprintf("%.6e", sum(values(I))), @sprintf("%.6e", exact),
                @sprintf("%.4f", maximum(r[h][:viol] for h in harm)))
    end
end
close(io)
println("✓ output/paper/bounds_results.csv")
