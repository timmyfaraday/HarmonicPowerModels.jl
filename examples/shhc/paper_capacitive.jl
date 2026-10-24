################################################################################
# HarmonicPowerModels.jl — SHHC paper §IV-C1: capacitive elements                #
#                                                                               #
# Networks: purely inductive base; a 19 Mvar capacitor bank at bus 29; the      #
# 33 kV lines as underground cables; the 132 kV lines as cables with 80 %       #
# reactor compensation (ieee30_scenarios.jl). One converter type, N = 21.        #
#   1. driving-point impedance |Z_30,30| over h = 1 … 15                         #
#   2. hosting capacity under arithmetic, robust, deterministic, exact 95 % and  #
#      IEC aggregation, and the equivalent summation exponent per harmonic       #
#   3. bank size at bus 29 from 0 to 20 Mvar                                     #
#                                                                               #
# Output: output/paper/cap_freqscan.csv, cap_results.csv, cap_size.csv          #
# Author: Tom Van Acker                                                          #
################################################################################

include(joinpath(@__DIR__, "paper_common.jl"))

datas = Dict(s => ieee30_scenario(s, ipopt) for s in PAPER_SCENARIOS)

io = open(joinpath(PAPER_OUT, "cap_freqscan.csv"), "w")
csvline(io, "h", PAPER_SCENARIOS...)
for h in 1.0:0.05:15.0
    zs = map(PAPER_SCENARIOS) do s
        Y, idx = ybus(datas[s], h); abs(inv(Y)[idx[30], idx[30]])
    end
    csvline(io, @sprintf("%.2f", h), [@sprintf("%.5e", z) for z in zs]...)
end
close(io)

io = open(joinpath(PAPER_OUT, "cap_results.csv"), "w")
csvline(io, "scenario", "arith", "rob", "det", "q95", "iec", "alpha5", "alpha7", "alpha11", "alpha13")
for s in PAPER_SCENARIOS
    r = recombination(s, NS, "local"; data = datas[s])
    csvline(io, s, [@sprintf("%.6e", tot(r, k)) for k in (:arith, :rob, :det, :q95, :iec)]..., alphas(r)...)
end
close(io)

io = open(joinpath(PAPER_OUT, "cap_size.csv"), "w")
csvline(io, "mvar", "q95", "rob", "det", "arith", "alpha5", "alpha7", "alpha11", "alpha13")
for mvar in (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5, 20.0)
    try
        d = mvar == 0 ? ieee30_scenario("inductive", ipopt) : ieee30_inductive_bank(PAPER_BANK_BUS, mvar / 100, ipopt)
        r = recombination("size$mvar", NS, "local"; data = d)
        csvline(io, mvar, [@sprintf("%.6e", tot(r, k)) for k in (:q95, :rob, :det, :arith)]..., alphas(r)...)
    catch err
        csvline(io, mvar, "NaN", "NaN", "NaN", "NaN", "NaN", "NaN", "NaN", "NaN")
    end
end
close(io)
println("✓ output/paper/cap_*.csv")
