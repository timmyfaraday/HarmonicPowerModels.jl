################################################################################
# HarmonicPowerModels.jl — Phase 2: recombination from mixed converter types,   #
# IEEE30, local angle reference                                                  #
#                                                                               #
# For every operable capacitive scenario, N = 8 and 21 harmonic units and the   #
# converter-type mixes of converter_types.jl, reports per harmonic the           #
# equivalent summation exponent of the 95% chance bound, and in total:           #
#   U       uplift of the 95% chance bound over the robust bound (I_95/I_rob)    #
#   gain    hosting with the deterministic angle over arithmetic summation       #
#   vs M0   95% hosting relative to the single-type population (same scenario)   #
#   IEC/95  hosting from the IEC 61000-3-6 summation law over the 95% value      #
# Binding constraints are the individual harmonic voltage limits (see            #
# ieee30_binding.txt); the same reduction as recomb_ieee30.jl.                   #
#                                                                               #
# Output: output/ieee30_phase2.txt                                              #
# Author: Tom Van Acker                                                          #
################################################################################

include(joinpath(@__DIR__, "recomb_ieee30.jl"))
include(joinpath(@__DIR__, "converter_types.jl"))

open(joinpath(@__DIR__, "output", "ieee30_phase2.txt"), "w") do io
    println(io, "IEEE30 Phase 2 — mixed converter types, local angle reference, binding: IHD.")
    println(io, "M0 all type A; M1 A/B(30°); M2 A/B(15°)/B(30°)/B(45°); M3 A/C (C = full-circle angle spread).")
    println(io, "α_eq: summation exponent reproducing the 95% hosting (IEC 1.4 for h=5,7; 2.0 for h=11,13); 'cap' = unit rating binds.\n")
    @printf(io, "%-10s %2s %-3s | %-23s | %7s %7s %7s %7s\n", "scenario", "N", "mix",
            "α_eq 95%  h=5/7/11/13", "U", "gain", "vs M0", "IEC/95")
    for scen in ["base", "banks", "cable33", "cable132c", "cableallc"], N in (8, 21)
        ref95 = NaN
        for mix in ["M0", "M1", "M2", "M3"]
            data  = ieee30_scenario(scen, ipopt)
            units = sort(parse.(Int, collect(keys(data["load"]))))[1:N]
            r = recombination(scen, N, "local"; ranges = mix_ranges(mix, units))
            tot(k) = sum(r[h][k] for h in harm)
            mix == "M0" && (ref95 = tot(:q95))
            @printf(io, "%-10s %2d %-3s | %s | %s %s %s %s\n", scen, N, mix,
                    join([fmt(r[h][:q95_alpha]) for h in harm], " "),
                    pct(tot(:q95), tot(:rob)), pct(tot(:det), tot(:arith)),
                    pct(tot(:q95), ref95), pct(tot(:iec), tot(:q95)))
            flush(io)
        end
    end
end
println("✓ ieee30_phase2.txt")
