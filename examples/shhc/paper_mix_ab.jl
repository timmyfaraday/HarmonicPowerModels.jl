################################################################################
# HarmonicPowerModels.jl — SHHC: mix of the diode (A) and thyristor (B) bridge   #
#                                                                               #
# Type A: §II, Table III. Type B: §II-D, Table IV, with the firing delay        #
# independent per unit (Table IV samples it jointly with L_ac and L_dc). The    #
# types alternate over the sorted units (11 A, 10 B for N = 21), each with its  #
# own angle ranges and Beta shape. Same networks and rules as                   #
# paper_capacitive.jl.                                                           #
#                                                                               #
# Output: output/paper_thyristor/mix_ab.csv (per scenario and harmonic)         #
# Author: Tom Van Acker                                                          #
################################################################################

include(joinpath(@__DIR__, "paper_common.jl"))

const SHAPE_A = (4.041, 4.776)
const SHAPE_B = (4.139, 3.659)
const RNG_A = Dict(5 => (4.75, 40.05), 7 => (-22.77, 27.51), 11 => (79.69, 164.01), 13 => (49.38, 147.67))
const RNG_B = Dict(5 => (-138.23, -10.87), 7 => (-265.07, -48.90), 11 => (-240.99, 44.65), 13 => (-337.52, 20.25))

units  = paper_units(ieee30_scenario("inductive", ipopt))
isB    = Dict(d => iseven(i) for (i, d) in enumerate(units))
shapes = Dict(d => isB[d] ? SHAPE_B : SHAPE_A for d in units)
ranges = Dict(h => Dict(d => deg2rad.(isB[d] ? RNG_B[h] : RNG_A[h]) for d in units) for h in harm)

out = joinpath(@__DIR__, "output", "paper_thyristor"); mkpath(out)
io = open(joinpath(out, "mix_ab.csv"), "w")
csvline(io, "scenario", "h", "arith", "rob", "det", "q95", "iec", "alpha_q95")
for s in PAPER_SCENARIOS
    r = recombination(s, NS, "local"; ranges = ranges, shapes = shapes, data = ieee30_scenario(s, ipopt))
    for h in harm
        csvline(io, s, h, [@sprintf("%.6e", r[h][k]) for k in (:arith, :rob, :det, :q95, :iec)]..., r[h][:q95_alpha])
    end
end
close(io)
println("✓ output/paper_thyristor/mix_ab.csv")
