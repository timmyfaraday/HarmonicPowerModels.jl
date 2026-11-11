################################################################################
# HarmonicPowerModels.jl — exact moments versus the expansion, harmonics 5–25   #
#                                                                               #
# sHHC_SOC with the Gaussian risk factor on the inductive IEEE 30-bus base      #
# network, N = 21 units, diode / thyristor / mix populations: squared-magnitude #
# moments from the expansion at δ = 2 and 3 ("moments" => "pce") and exact      #
# ("moments" => "exact"), against the exact 95 % quantile (Monte Carlo), with   #
# the violation probability of the binding individual harmonic voltage limit.   #
#                                                                               #
# Output: output/paper_h25/exact_moments.csv                                    #
# Author: Tom Van Acker                                                          #
################################################################################

ENV["SHHC_HMAX"] = "25"
include(joinpath(@__DIR__, "paper_common.jl"))
data  = ieee30_scenario("inductive", ipopt)
units = paper_units(data)
out = open(joinpath(PAPER_OUT, "exact_moments.csv"), "w"); println(out, "population,method,h,I_soc,I_mc,viol,time")
for pop in ("diode", "thyristor", "mix")
    ranges, shapes = population(units, pop)
    ref = recombination("inductive", NS, "local"; data = data, ranges = ranges, shapes = shapes)
    for (meth, δ, mom) in (("pce2", 2, "pce"), ("pce3", 3, "pce"), ("exact", 2, "exact"))
        p = paper_pce(units, nothing; bound = "gaussian", pop = pop); p["delta"] = δ; p["moments"] = mom
        t = @elapsed sol = HPM.solve_shhc(restricted(data, units), sHHC_SOC, mosek; pce_params = p, hpf_optimizer = ipopt)
        I = Dict(h => sol["solution"][h] for h in harm)
        v = recombination("inductive", NS, "local"; data = data, ranges = ranges, shapes = shapes, eval_I = I)
        for h in harm; println(out, join((pop, meth, h, I[h], ref[h][:q95], v[h][:viol], t), ",")); flush(out); end
        println("done $pop $meth $(sol["termination_status"]) $(round(t, digits = 2)) s")
    end
end
close(out)
