################################################################################
# HarmonicPowerModels.jl — SHHC paper §IV-B: computational tractability          #
#                                                                               #
# Wall time of sHHC_SOC with the empirical risk factor (fundamental power      #
# flow, polynomial chaos data, the sampling calibration of λ* and every cone    #
# program of the calibration loop) and the solver time of Mosek alone, summed   #
# over the loop, on the IEEE 30-bus networks of the paper, against             #
#   1. the number of harmonic units N = 4, 8, 13, 21 (diode bridges);           #
#   2. the converter population (N = 21);                                       #
#   3. the number of harmonics |H| = 1, 2, 3, 4 (N = 21, diode bridges);        #
#   4. the two networks of §IV-B1 (N = 21, diode bridges).                      #
# Each case is solved twice with an empty basis cache and the second run is    #
# reported, so that Julia's compilation is excluded; a third run reuses the    #
# cached basis (wall_cached_s). The cardinality |K| of the basis is reported.   #
#                                                                               #
# Output: output/paper/tractability.csv                                         #
# Author: Tom Van Acker                                                          #
################################################################################

include(joinpath(@__DIR__, "paper_common.jl"))

function timed_soc(data, units, pop, hs)
    hd() = begin
        d = HPM.replicate(deepcopy(data), H = [1; hs]); d["angle_reference"] = "local"
        for h in hs, k in collect(keys(d["nw"]["$h"]["load"]))
            parse(Int, k) in units || delete!(d["nw"]["$h"]["load"], k)
        end
        d
    end
    p = paper_pce(units, nothing; bound = "gaussian", pop = pop)   # gaussian: first solve of the loop
    p["harmonics"] = Dict{Int,Any}(h => p["harmonics"][h] for h in hs)
    solve(d) = HPM.solve_shhc_calibrated(d, sHHC_SOC, mosek; pce_params = deepcopy(p), hpf_optimizer = ipopt,
                                         nsamp = 20000)
    local sol, t
    for _ in 1:2                       # compile, then a cold run (basis built)
        empty!(HPM._PCEBASIS_CACHE)
        t = @elapsed sol = solve(hd())
    end
    tc = @elapsed solve(hd())          # warm run: basis taken from the cache
    # solver time of one cone program, times the number of solves in the loop
    d1 = hd(); s1 = HPM.solve_shhc(d1, sHHC_SOC, mosek; pce_params = deepcopy(p), hpf_optimizer = ipopt)
    nsolve = sol["iterations"] + 1
    return t, nsolve * s1["solve_time"], d1["pce"]["basis"].dim, string(sol["termination_status"]), tc, nsolve
end

io = open(joinpath(PAPER_OUT, "tractability.csv"), "w")
csvline(io, "study", "network", "buses", "N", "population", "H", "K", "wall_s", "solver_s", "status", "wall_cached_s", "solves")
row(study, net, data, N, pop, hs, r) = (csvline(io, study, net, length(data["bus"]), N, pop, length(hs), r[3],
                                               @sprintf("%.3f", r[1]), @sprintf("%.4f", r[2]), r[4], @sprintf("%.3f", r[5]), r[6]); flush(io))
base = ieee30_scenario("inductive", ipopt)
for N in (4, 8, 13, 21)
    u = paper_units(base, N); row("units", "base", base, N, "diode", harm, timed_soc(base, u, "diode", harm))
end
u21 = paper_units(base)
for pop in POPULATIONS
    row("populations", "base", base, NS, pop, harm, timed_soc(base, u21, pop, harm))
end
for n in 1:4
    hs = harm[1:n]; row("harmonics", "base", base, NS, "diode", hs, timed_soc(base, u21, "diode", hs))
end
for s in ("inductive", "inductive_cable33")
    d = ieee30_scenario(s, ipopt); row("networks", s, d, NS, "diode", harm, timed_soc(d, paper_units(d), "diode", harm))
end
close(io)
println("✓ output/paper/tractability.csv")
