################################################################################
# HarmonicPowerModels.jl — SHHC paper §IV-C2: converter types                    #
#                                                                               #
# T = 1, 2, 4, 8, 21 converter types among N = 21 units, mean firing delays      #
# spread evenly over [0°, 60°), assigned cyclically; type 1 is the §II converter #
# and every other type has one random variable for its firing delay, shared by  #
# its units (paper_common.jl). T = 21 gives every unit its own type. Networks:  #
# purely inductive base and the base with the capacitor bank at bus 29.          #
# Reported: hosting capacity under arithmetic, robust, deterministic, exact 95 % #
# and IEC aggregation, and the equivalent summation exponent per harmonic.       #
#                                                                               #
# Output: output/paper/types_results.csv                                         #
# Author: Tom Van Acker                                                          #
################################################################################

include(joinpath(@__DIR__, "paper_common.jl"))

io = open(joinpath(PAPER_OUT, "types_results.csv"), "w")
csvline(io, "scenario", "T", "arith", "rob", "det", "q95", "iec", "alpha5", "alpha7", "alpha11", "alpha13")
for s in ("inductive", "inductive_bank")
    data  = ieee30_scenario(s, ipopt)
    units = paper_units(data)
    for T in (1, 2, 4, 8, 21)
        r = recombination(s, NS, "local"; data = data, types = types_spread(units, T))
        csvline(io, s, T, [@sprintf("%.6e", tot(r, k)) for k in (:arith, :rob, :det, :q95, :iec)]..., alphas(r)...)
    end
end
close(io)
println("✓ output/paper/types_results.csv")
