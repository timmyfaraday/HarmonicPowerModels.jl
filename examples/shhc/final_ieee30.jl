################################################################################
# HarmonicPowerModels.jl — final IEEE30 results for the SHHC paper               #
#                                                                               #
# Scenarios base, banks, cable33, cable132c; N = 8 and 21 harmonic units;        #
# local angle reference; two populations:                                        #
#   A   all units type A (the §II converter)                                     #
#   AB  50/50 type A and type B, type B with a firing delay of 30° mean and 30°  #
#       width shared by all type-B units (one type germ — working assumption)    #
# Per case: the reduced hosting (binding IHD constraints) under arithmetic      #
# summation, robust, deterministic, exact 95% and the IEC law, the equivalent   #
# summation exponents, and the full sHHC_SOC hosting for the Gaussian, Cantelli #
# and empirical bounds (Mosek), plus dHHC at E[φ] as a check on the reduction.  #
#                                                                               #
# Output: output/ieee30_final.csv, output/ieee30_final.txt                      #
# Author: Tom Van Acker                                                          #
################################################################################

include(joinpath(@__DIR__, "recomb_ieee30.jl"))
include(joinpath(@__DIR__, "converter_types.jl"))

function restricted(data, units)
    hd = HPM.replicate(deepcopy(data), H = H); hd["angle_reference"] = "local"
    for h in harm, k in collect(keys(hd["nw"]["$h"]["load"]))
        parse(Int, k) in units || delete!(hd["nw"]["$h"]["load"], k)
    end
    return hd
end

csv = open(joinpath(@__DIR__, "output", "ieee30_final.csv"), "w")
txt = open(joinpath(@__DIR__, "output", "ieee30_final.txt"), "w")
println(csv, "scenario,N,population,arith,rob,det,q95,iec,soc_gauss,soc_cant,soc_emp,dhhc,alpha5,alpha7,alpha11,alpha13")
println(txt, "IEEE30 final results — local angle reference, Mosek; Σ_h I_h [p.u.]; AB = 50/50 A and B(30° ± 15°, shared germ).\n")
@printf(txt, "%-10s %2s %-2s | %9s %9s %9s %9s %9s | %9s %9s %9s | %-23s\n", "scenario", "N", "", "arith", "robust",
        "det", "q95", "IEC", "SOC gau", "SOC cant", "SOC emp", "α_eq 95% h=5/7/11/13")
for scen in ["base", "banks", "cable33", "cable132c"], N in (8, 21), pop in ("A", "AB")
    data  = ieee30_scenario(scen, ipopt)
    units = sort(parse.(Int, collect(keys(data["load"]))))[1:N]
    types = pop == "AB" ? typeB_germ(units) : nothing
    r = recombination(scen, N, "local"; types = types)
    tot(k) = sum(r[h][k] for h in harm)
    # full stochastic model
    p0 = Dict{String,Any}(deepcopy(PCE_PARAMS)...)
    p0["harmonics"] = Dict{Int,Any}(h => rng[h] for h in harm); p0["harmonic_units"] = units
    types === nothing || (p0["type_germs"] = types)
    soc = Dict{String,Float64}()
    for b in ("gaussian", "cantelli")
        p = deepcopy(p0); p["bound"] = b
        s = HPM.solve_shhc(restricted(data, units), sHHC_SOC, mosek; pce_params = p, hpf_optimizer = ipopt)
        soc[b] = sum(s["solution"][h] for h in harm)
    end
    se = HPM.solve_shhc_calibrated(restricted(data, units), sHHC_SOC, mosek;
                                   pce_params = deepcopy(p0), hpf_optimizer = ipopt, nsamp = 20000)
    soc["empirical"] = sum(se["solution"][h] for h in harm)
    # dHHC at E[φ] (type germ at its mean) — must equal the reduced deterministic value
    hd = restricted(data, units)
    for h in harm, d in units
        lo, hi = rng[h]; φ = lo + (hi - lo) * Eω
        if types !== nothing && d in types["B"]["units"]
            sl, sh = types["B"]["shift"][h]; φ += sl + (sh - sl) * Eω
        end
        hd["nw"]["$h"]["load"]["$d"]["ref_angle"] = rad2deg(φ)
    end
    sd = HPM.solve_hhc(hd, dHHC_SOC, mosek, ipopt)
    dh = sum(sd["solution"]["nw"]["$h"]["load"]["$(units[1])"]["cmd"] for h in harm)
    a = [r[h][:q95_alpha] for h in harm]
    @printf(csv, "%s,%d,%s,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%s\n", scen, N, pop,
            tot(:arith), tot(:rob), tot(:det), tot(:q95), tot(:iec), soc["gaussian"], soc["cantelli"],
            soc["empirical"], dh, join(string.(a), ","))
    @printf(txt, "%-10s %2d %-2s | %.3e %.3e %.3e %.3e %.3e | %.3e %.3e %.3e | %s   (dHHC %+.1e%%)\n", scen, N, pop,
            tot(:arith), tot(:rob), tot(:det), tot(:q95), tot(:iec), soc["gaussian"], soc["cantelli"], soc["empirical"],
            join([fmt(x) for x in a], " "), 100(dh / tot(:det) - 1))
    flush(csv); flush(txt)
end
close(csv); close(txt)
println("✓ ieee30_final.{csv,txt}")
