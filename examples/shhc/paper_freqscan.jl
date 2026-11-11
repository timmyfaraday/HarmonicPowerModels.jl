################################################################################
# HarmonicPowerModels.jl — SHHC paper §IV: frequency scan up to the 25th        #
#                                                                               #
# Driving-point impedance |Z_ii| at bus 30 (end of the longest 33 kV feeder)    #
# and bus 26 (binding at the 19th to 25th harmonic with the 33 kV cables) for   #
# the inductive base network and with all 33 kV lines as cables, h = 1 … 26.   #
#                                                                               #
# Output: output/paper[_h25]/freqscan.csv                                       #
# Author: Tom Van Acker                                                          #
################################################################################

include(joinpath(@__DIR__, "paper_common.jl"))

nets  = ["inductive", "inductive_cable33"]
datas = Dict(s => ieee30_scenario(s, ipopt) for s in nets)
io = open(joinpath(PAPER_OUT, "freqscan.csv"), "w")
csvline(io, "h", "base30", "cable30", "base26", "cable26")
for h in 1.0:0.05:26.0
    z = Dict((s, b) => (Y = ybus(datas[s], h); abs(inv(Y[1])[Y[2][b], Y[2][b]])) for s in nets, b in (30, 26))
    csvline(io, @sprintf("%.2f", h), [@sprintf("%.5e", z[(s, b)]) for b in (30, 26) for s in nets]...)
end
close(io)
println("✓ ", joinpath(PAPER_OUT, "freqscan.csv"))
