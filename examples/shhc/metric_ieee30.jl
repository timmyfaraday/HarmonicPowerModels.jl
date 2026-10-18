################################################################################
# HarmonicPowerModels.jl — reference metric for the stochastic uplift, IEEE30    #
#                                                                               #
# Per harmonic, with the IHD constraints binding (see ieee30_binding.txt), the  #
# hosting current under absolute equality is the tightest bus:                   #
#   I_det,h = min_i  L_ih / |Σ_j c_ij(E[φ])|        deterministic, E[φ]          #
#   I_rob,h = min_i  L_ih / max_φ∈support |Σ_j c_ij| robust over the support     #
#   I_95,h  = min_i  L_ih / q95(|Σ_j c_ij|)          exact per-constraint 95%    #
# with L_ih = ihdmax·|V_i,1| and c_ij = Z_h[i, bus_j] e^{jφ_j}. The unit rating   #
# caps every value. Reported totals Σ_h and the uplift U = I_95 / I_rob.        #
#                                                                               #
# Output: output/ieee30_metric.txt                                              #
# Author: Tom Van Acker                                                          #
################################################################################

include(joinpath(@__DIR__, "diag_binding_ieee30.jl"))

function metric(scen, N, aref)
    data = ieee30_scenario(scen, ipopt)
    units = sort(parse.(Int, collect(keys(data["load"]))))[1:N]
    hd = HPM.replicate(deepcopy(data), H = H)
    hd["angle_reference"] = aref
    for h in harm, k in collect(keys(hd["nw"]["$h"]["load"]))
        parse(Int, k) in units || delete!(hd["nw"]["$h"]["load"], k)
    end
    HPM.update_hdata_with_fundamental_hpf_results!(hd, dHHC_NLP, ipopt)
    Z, bus_pos = HPM.reduced_harmonic_impedance(hd)
    f = hd["nw"]["1"]
    ubus = Dict(d => f["load"]["$d"]["load_bus"] for d in units)
    beta = Beta(α̂, β̂); Random.seed!(1)
    tot = Dict(:det => 0.0, :rob => 0.0, :q95 => 0.0)
    for h in harm
        off = [deg2rad(get(hd["nw"]["$h"]["load"]["$d"], "angle_offset", 0.0)) for d in units]
        lo, hi = rng[h]; lo_j = lo .+ off; hi_j = hi .+ off
        φE = lo_j .+ (hi - lo) * Eω
        W  = [lo_j .+ (hi - lo) .* rand(beta, N) for _ in 1:20_000]     # shared angle samples
        cap = minimum(hd["nw"]["$h"]["load"]["$d"]["c_rating"] for d in units)
        Id = Ir = Iq = cap
        for (b, bus) in f["bus"]
            i = parse(Int, b); haskey(bus_pos, i) || continue
            L = hd["nw"]["$h"]["bus"][b]["ihdmax"] * bus["vm"]
            z = [Z[h][bus_pos[i], bus_pos[ubus[d]]] for d in units]
            Id = min(Id, L / abs(sum(z .* cis.(φE))))
            Ir = min(Ir, L / robust_abs(z, lo_j, hi_j))
            Iq = min(Iq, L / quantile([abs(sum(z .* cis.(w))) for w in W], 0.95))
        end
        tot[:det] += Id; tot[:rob] += Ir; tot[:q95] += Iq
    end
    return tot
end

open(joinpath(@__DIR__, "output", "ieee30_metric.txt"), "w") do io
    println(io, "IEEE30 — IHD-based hosting Σ_h I_h [p.u.]: deterministic at E[φ], robust over the angle")
    println(io, "support, exact per-constraint 95% quantile; U = I_95/I_rob, D = I_95/I_det.\n")
    @printf(io, "%-10s %3s %-6s | %10s %10s %10s | %7s %7s\n", "scenario", "N", "ref", "I_det", "I_rob", "I_95", "U", "D")
    for scen in ["base", "banks", "cable132c"], N in (8, 21), aref in ("global", "local")
        t = metric(scen, N, aref)
        @printf(io, "%-10s %3d %-6s | %10.4e %10.4e %10.4e | %+6.1f%% %+6.1f%%\n", scen, N, aref,
                t[:det], t[:rob], t[:q95], 100(t[:q95]/t[:rob]-1), 100(t[:q95]/t[:det]-1))
        flush(io)
    end
end
println("✓ ieee30_metric.txt")
