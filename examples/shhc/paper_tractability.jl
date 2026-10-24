################################################################################
# HarmonicPowerModels.jl — SHHC paper §IV-B: computational tractability          #
#                                                                               #
# Wall time of sHHC_SOC (fundamental power flow, polynomial chaos data,          #
# precomputation of the coefficients α and the cone program) and the solver     #
# time of Mosek alone, on the IEEE 30-bus networks of the paper, against        #
#   1. the number of harmonic units N = 4, 8, 13, 21 (one converter type);      #
#   2. the number of converter types T = 1, 2, 4, 8, 21 (N = 21);               #
#   3. the number of harmonics |H| = 1, 2, 3, 4 (N = 21, one type);             #
#   4. the four networks of §IV-B1 (N = 21, one type).                          #
# Each case is solved twice with an empty basis cache and the second run is    #
# reported, so that Julia's compilation is excluded; a third run reuses the    #
# cached basis (wall_cached_s). The cardinality |K| of the basis is reported.   #
#                                                                               #
# Output: output/paper/tractability.csv                                         #
# Author: Tom Van Acker                                                          #
################################################################################

include(joinpath(@__DIR__, "paper_common.jl"))

function timed_soc(data, units, types, hs)
    hd() = begin
        d = HPM.replicate(deepcopy(data), H = [1; hs]); d["angle_reference"] = "local"
        for h in hs, k in collect(keys(d["nw"]["$h"]["load"]))
            parse(Int, k) in units || delete!(d["nw"]["$h"]["load"], k)
        end
        d
    end
    p = paper_pce(units, types; bound = "gaussian")
    p["harmonics"] = Dict{Int,Any}(h => rng[h] for h in hs)
    local sol, t, d
    for _ in 1:2                       # compile, then a cold run (basis built)
        empty!(HPM._PCEBASIS_CACHE); d = hd()
        t = @elapsed sol = HPM.solve_shhc(d, sHHC_SOC, mosek; pce_params = deepcopy(p), hpf_optimizer = ipopt)
    end
    d2 = hd()                          # warm run: basis taken from the cache
    tc = @elapsed HPM.solve_shhc(d2, sHHC_SOC, mosek; pce_params = deepcopy(p), hpf_optimizer = ipopt)
    return t, sol["solve_time"], d["pce"]["basis"].dim, string(sol["termination_status"]), tc
end

io = open(joinpath(PAPER_OUT, "tractability.csv"), "w")
csvline(io, "study", "network", "buses", "N", "T", "H", "K", "wall_s", "solver_s", "status", "wall_cached_s")
row(study, net, data, N, T, hs, r) = (csvline(io, study, net, length(data["bus"]), N, T, length(hs), r[3],
                                             @sprintf("%.3f", r[1]), @sprintf("%.4f", r[2]), r[4], @sprintf("%.3f", r[5])); flush(io))
base = ieee30_scenario("inductive", ipopt)
for N in (4, 8, 13, 21)
    u = paper_units(base, N); row("units", "base", base, N, 1, harm, timed_soc(base, u, nothing, harm))
end
u21 = paper_units(base)
for T in (1, 2, 4, 8, 21)
    row("types", "base", base, NS, T, harm, timed_soc(base, u21, types_spread(u21, T), harm))
end
for n in 1:4
    hs = harm[1:n]; row("harmonics", "base", base, NS, 1, hs, timed_soc(base, u21, nothing, hs))
end
for s in PAPER_SCENARIOS
    d = ieee30_scenario(s, ipopt); row("networks", s, d, NS, 1, harm, timed_soc(d, paper_units(d), nothing, harm))
end
close(io)
println("✓ output/paper/tractability.csv")
