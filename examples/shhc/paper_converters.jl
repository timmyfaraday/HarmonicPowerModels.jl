################################################################################
# HarmonicPowerModels.jl — SHHC paper §IV-B: converter populations               #
#                                                                               #
# Networks: the purely inductive base and the 33 kV lines as underground cables #
# (ieee30_scenarios.jl). Populations over the N = 21 units: all diode bridges   #
# (Table III), all thyristor bridges (Table IV, firing delay independent per    #
# unit) and a mix of both (alternating). Hosting capacity under arithmetic,     #
# robust, mean-angle, exact 95 % and IEC aggregation, and the equivalent        #
# summation exponent per harmonic.                                               #
#                                                                               #
# Output: output/paper/converters_results.csv                                   #
# Author: Tom Van Acker                                                          #
################################################################################

include(joinpath(@__DIR__, "paper_common.jl"))

const NETWORKS = ["inductive", "inductive_cable33"]

io = open(joinpath(PAPER_OUT, "converters_results.csv"), "w")
csvline(io, "network", "population", "arith", "rob", "det", "q95", "iec", ["alpha$h" for h in harm]...)
for s in NETWORKS
    data  = ieee30_scenario(s, ipopt)
    units = paper_units(data)
    for pop in POPULATIONS
        ranges, shapes = population(units, pop)
        r = recombination(s, NS, "local"; data = data, ranges = ranges, shapes = shapes)
        csvline(io, s, pop, [@sprintf("%.6e", tot(r, k)) for k in (:arith, :rob, :det, :q95, :iec)]..., alphas(r)...)
    end
end
close(io)
println("✓ ", joinpath(PAPER_OUT, "converters_results.csv"))
