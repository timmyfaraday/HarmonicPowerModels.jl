################################################################################
# HarmonicPowerModels.jl — Phase 2 sensitivity sweep, IEEE30, local reference   #
#                                                                               #
# On base and banks at N = 8 and 21:                                            #
#   1. firing delay αf = 0 … 60° for a 50/50 mix of A and B(αf)                  #
#   2. share of type B(30°) from 0 to 100%                                      #
#   3. firing delay as a range: B(30°, w) with w = 0 … 60°, 50/50 mix           #
# Reported: equivalent summation exponent of the 95% hosting per harmonic, the  #
# uplift over the robust bound U, and the 95% hosting relative to all type A.   #
#                                                                               #
# Output: output/ieee30_sweep.txt and output/ieee30_sweep.csv                   #
# Author: Tom Van Acker                                                          #
################################################################################

include(joinpath(@__DIR__, "recomb_ieee30.jl"))
include(joinpath(@__DIR__, "converter_types.jl"))

txt = open(joinpath(@__DIR__, "output", "ieee30_sweep.txt"), "w")
csv = open(joinpath(@__DIR__, "output", "ieee30_sweep.csv"), "w")
println(csv, "sweep,scenario,N,value,alpha5,alpha7,alpha11,alpha13,U,vs_A,q95,rob,det,arith")
println(txt, "IEEE30 Phase 2 sensitivity sweep — local angle reference, binding: IHD.")
println(txt, "α_eq: summation exponent reproducing the 95% hosting (IEC 1.4 for h=5,7; 2.0 for h=11,13).")
println(txt, "U: 95% hosting over the robust bound; vs A: 95% hosting over the all-type-A population.\n")

for scen in ["base", "banks"], N in (8, 21)
    data  = ieee30_scenario(scen, ipopt)
    units = sort(parse.(Int, collect(keys(data["load"]))))[1:N]
    refA  = recombination(scen, N, "local"; ranges = share_ranges(0.0, :A, units))
    q95A  = sum(refA[h][:q95] for h in harm)
    function row(sweep, label, rg)
        r = recombination(scen, N, "local"; ranges = rg)
        tot(k) = sum(r[h][k] for h in harm)
        @printf(txt, "  %-12s | %s | %8s %8s\n", label, join([fmt(r[h][:q95_alpha]) for h in harm], " "),
                pct(tot(:q95), tot(:rob)), pct(tot(:q95), q95A))
        @printf(csv, "%s,%s,%d,%s,%s,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n", sweep, scen, N, label,
                join([string(r[h][:q95_alpha]) for h in harm], ","),
                tot(:q95) / tot(:rob), tot(:q95) / q95A, tot(:q95), tot(:rob), tot(:det), tot(:arith))
        flush(txt); flush(csv)
    end
    println(txt, repeat('═', 70)); println(txt, "$scen  N = $N"); println(txt, repeat('═', 70))
    @printf(txt, "  %-12s | %-23s | %8s %8s\n", "", "α_eq  h=5/7/11/13", "U", "vs A")
    println(txt, "  1. firing delay, 50/50 A and B(αf)")
    for αf in 0:5:60
        row("delay", "αf=$(αf)°", share_ranges(0.5, (:B, Float64(αf)), units))
    end
    println(txt, "  2. share of B(30°)")
    for s in (0.0, 0.25, 0.5, 0.75, 1.0)
        row("share", "share=$(round(Int, 100s))%", share_ranges(s, (:B, 30.0), units))
    end
    println(txt, "  3. firing delay range B(30°, w), 50/50")
    for w in (0, 10, 20, 30, 40, 60)
        row("range", "w=$(w)°", share_ranges(0.5, (:B, 30.0, Float64(w)), units))
    end
    println(txt)
end
close(txt); close(csv)
println("✓ ieee30_sweep.{txt,csv}")
