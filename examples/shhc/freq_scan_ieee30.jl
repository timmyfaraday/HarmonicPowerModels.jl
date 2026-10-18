################################################################################
# HarmonicPowerModels.jl — frequency scan of the IEEE30 capacitive scenarios    #
#                                                                               #
# For every scenario in ieee30_scenarios.jl:                                    #
#   1. fundamental power flow — voltage range and maximum branch loading, which #
#      decide whether the scenario can be operated at all;                      #
#   2. driving-point impedance Z_ii(h) on a fine harmonic grid h = 1 … 25, using #
#      the same stamps as the sHHC_SOC admittance (branch π-model with tap,      #
#      bus shunts, generator subtransient admittance), checked against           #
#      `_harmonic_admittance` at the characteristic harmonics;                   #
#   3. at h = 5, 7, 11, 13: the spread of the transfer-impedance angles arg Z_ij #
#      from the harmonic units j to each bus i, i.e., how far the network alone  #
#      already rotates the contributions apart.                                  #
#                                                                               #
# Output: output/ieee30_freqscan.csv (scenario,bus,h,absZ,angZ) and             #
#         output/ieee30_freqscan.txt (summary).                                  #
# Author: Tom Van Acker                                                          #
################################################################################

using HarmonicPowerModels, PowerModels, Ipopt, LinearAlgebra, Printf
const HPM = HarmonicPowerModels
PowerModels.silence()
include(joinpath(@__DIR__, "ieee30_scenarios.jl"))

ipopt  = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
outdir = joinpath(@__DIR__, "output"); mkpath(outdir)
hgrid  = 1.0:0.05:25.0
harm   = [5, 7, 11, 13]
scan_buses = [4, 10, 17, 24, 30]

"Harmonic admittance at a (possibly non-integer) harmonic order h from fundamental data."
function ybus(data, h::Real)
    ids, idx = HPM._bus_id_map(data)
    Y = zeros(ComplexF64, length(ids), length(ids))
    for br in values(data["branch"])
        get(br, "br_status", 1) == 0 && continue
        f, t = idx[br["f_bus"]], idx[br["t_bus"]]
        tm = get(br, "tap", 1.0); tm == 0.0 && (tm = 1.0)
        tc = tm * cis(get(br, "shift", 0.0))
        ys = 1 / (br["br_r"] * sqrt(h) + im * br["br_x"] * h)
        yf = get(br, "g_fr", 0.0) + im * get(br, "b_fr", 0.0) * h
        yt = get(br, "g_to", 0.0) + im * get(br, "b_to", 0.0) * h
        Y[f, f] += (ys + yf) / abs2(tc); Y[t, t] += ys + yt
        Y[f, t] -= ys / conj(tc);        Y[t, f] -= ys / tc
    end
    for sh in values(get(data, "shunt", Dict{String,Any}()))
        bs = sh["bs"]
        Y[idx[sh["shunt_bus"]], idx[sh["shunt_bus"]]] += sh["gs"] + im * (bs >= 0 ? bs * h : bs / h)
    end
    for g in values(data["gen"])
        rsc, xsc = get(g, "rsc", 0.0), get(g, "xsc", 0.0)
        (iszero(rsc) && iszero(xsc)) && continue
        Y[idx[g["gen_bus"]], idx[g["gen_bus"]]] += 1 / (rsc * sqrt(h) + im * xsc * h)
    end
    return Y, idx
end

wrapdeg(a) = mod(a + 180.0, 360.0) - 180.0
"Smallest arc containing all angles [deg]."
function angle_spread(angs)
    a = sort(mod.(angs, 360.0))
    gaps = vcat(diff(a), a[1] + 360.0 - a[end])
    return 360.0 - maximum(gaps)
end

csv = open(joinpath(outdir, "ieee30_freqscan.csv"), "w")
txt = open(joinpath(outdir, "ieee30_freqscan.txt"), "w")
println(csv, "scenario,bus,h,absZ,angZ")
sep() = println(txt, repeat('═', 92))

for name in ["base", "banks", "cable33", "cable132", "cableall", "cable132c", "cableallc"]
    data = try
        ieee30_scenario(name, ipopt)
    catch err
        sep(); println(txt, "SCENARIO $(uppercase(name))"); sep()
        println(txt, "Not operable: no AC-OPF dispatch respects the branch ratings ($(sprint(showerror, err)))\n")
        continue
    end
    unit_buses = Dict(n => [data["load"]["$l"]["load_bus"] for l in sort(parse.(Int, collect(keys(data["load"]))))[1:n]] for n in (8, 21))

    # 1. fundamental power flow
    hd = HPM.replicate(deepcopy(data), H = [1, 5, 7, 11, 13])
    pf_ok = true
    vmin = vmax = load_max = NaN; worst = ""
    try
        HPM.update_hdata_with_fundamental_hpf_results!(hd, dHHC_NLP, ipopt)
        f   = hd["nw"]["1"]
        vm  = [b["vm"] for b in values(f["bus"])]
        vmin, vmax = extrema(vm)
        load_max, k = findmax(Dict(k => max(b["cm_fr"], b["cm_to"]) / b["c_rating"] for (k, b) in f["branch"]))
        worst = "$(f["branch"][k]["f_bus"])-$(f["branch"][k]["t_bus"])"
    catch err
        pf_ok = false
    end

    # 2. consistency check against the model's own admittance at the harmonics
    _, bus_idx = HPM._bus_id_map(data)
    dev = maximum(maximum(abs.(ybus(data, h)[1] .- HPM._harmonic_admittance(hd["nw"]["$h"], bus_idx))) for h in harm)

    # driving-point scan
    peaks = Dict{Int,Vector{Float64}}()
    for b in scan_buses
        zs = Float64[]
        for h in hgrid
            Y, idx = ybus(data, h)
            z = inv(Y)[idx[b], idx[b]]
            push!(zs, abs(z))
            @printf(csv, "%s,%d,%.2f,%.6e,%.3f\n", name, b, h, abs(z), rad2deg(angle(z)))
        end
        peaks[b] = [hgrid[i] for i in 2:length(zs)-1 if zs[i] > zs[i-1] && zs[i] > zs[i+1]]
    end

    # 3. transfer-impedance angle spread at the characteristic harmonics
    sep(); println(txt, "SCENARIO $(uppercase(name))"); sep()
    if pf_ok
        @printf(txt, "Fundamental: |V| in [%.4f, %.4f] p.u., max branch loading %.1f%% (branch %s)%s\n",
                vmin, vmax, 100load_max, worst, load_max >= 1 ? "  ← OVERLOADED: no harmonic headroom" : "")
        @printf(txt, "Dispatch: %s\n", join([@sprintf("gen@%d Pg=%.3f Vg=%.4f", g["gen_bus"], g["pg"], g["vg"])
                                            for (_, g) in sort(collect(data["gen"]), by = x -> parse(Int, x[1]))], ", "))
    else
        println(txt, "Fundamental: power flow did not solve")
    end
    @printf(txt, "Admittance check vs _harmonic_admittance at h=5,7,11,13: max |ΔY| = %.2e\n", dev)
    startswith(name, "cable") &&
        @printf(txt, "Cable length: %.0f km in total\n", sum(get(b, "cable_km", 0.0) for b in values(data["branch"])))
    println(txt, "\nParallel-resonance peaks of |Z_ii| (harmonic order):")
    for b in scan_buses
        @printf(txt, "  bus %-3d %s\n", b, join([@sprintf("%.2f", p) for p in peaks[b]], ", "))
    end
    for n in (8, 21)
        println(txt, "\nTransfer-impedance angle spread over the $n units' buses [deg]  (bus with the largest Σ_j|Z_ij| / max over buses):")
        for h in harm
            Y, idx = ybus(data, h); Z = inv(Y)
            cols = [idx[b] for b in unit_buses[n]]
            exposure = [sum(abs.(Z[i, cols])) for i in 1:size(Z, 1)]
            iexp = argmax(exposure)
            spreads = [angle_spread(rad2deg.(angle.(Z[i, cols]))) for i in 1:size(Z, 1)]
            inv_idx = Dict(v => k for (k, v) in idx)
            @printf(txt, "  h=%-3d exposed bus %-3d spread %6.1f   |  max spread %6.1f at bus %d\n",
                    h, inv_idx[iexp], spreads[iexp], maximum(spreads), inv_idx[argmax(spreads)])
        end
    end
    println(txt)
end
close(csv); close(txt)
println("✓ ieee30_freqscan.{csv,txt}")
