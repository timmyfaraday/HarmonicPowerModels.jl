################################################################################
# Diagnostic: for the constructive baseline, is the IMPEDANCE angle shift        #
# negligible vs the transformer VECTOR-GROUP rotation?                           #
#                                                                               #
# For the binding bus n* and each harmonic unit d, the angle that aligns d's     #
# current contribution constructively at n* is  φ_d = const − arg(Z_{n*,d}).     #
# arg(Z_{n*,d}) = [vector-group rotation] + [impedance angle].                   #
# We report both parts so the user can judge whether ignoring impedance is OK.   #
################################################################################
using HarmonicPowerModels, PowerModels
using LinearAlgebra, Printf
const HPM = HarmonicPowerModels
const PMs = PowerModels

include(joinpath(@__DIR__, "pce_params.jl"))
data = PMs.parse_file(joinpath(@__DIR__, "data", "industrial_network_shhc.m"))
H = [1,5,7,11,13]; harm = [5,7,11,13]
hd = HPM.replicate(deepcopy(data), H = H)

fund = hd["nw"]["1"]
bus_ids = sort([b["bus_i"] for (_,b) in fund["bus"]])
bus_idx = Dict(id=>i for (i,id) in enumerate(bus_ids))
Nn = length(bus_ids)
ref_ids = Set(bus_idx[b["bus_i"]] for (_,b) in fund["bus"] if b["bus_type"]==3)
non_ref = [i for i in 1:Nn if i ∉ ref_ids]
pos = Dict(id=>k for (k,id) in enumerate(non_ref))

# unit buses (loads with nonzero Pd in fundamental)
unit_buses = sort([l["load_bus"] for (_,l) in fund["load"]])
# vector-group cumulative rotation per bus (from bus_harmonics ref_angle), deg
vg_ang = Dict{Int,Float64}()
for (_,b) in fund["bus"]
    vg_ang[b["bus_i"]] = get(b, "ref_angle", 0.0)
end
seqsign(h) = (h % 3 == 1) ? +1.0 : (h % 3 == 2) ? -1.0 : 0.0   # pos:+, neg:−, zero:0

nstar = bus_idx[8]   # binding bus from the recombination diagnostic (THD@bus8)

println("Binding bus n* = 8.  Per harmonic, per unit d:")
println("  argZ = arg(Z_{8,d}) [full, deg];  VG = seq·(ref_angle_d − ref_angle_8);")
println("  imp  = argZ − VG  (the impedance angle shift you propose to ignore).\n")
@printf("%-4s %-8s %-10s %-10s %-10s %-10s\n","h","unit@bus","|Z|","argZ[°]","VG[°]","imp[°]")
for h in harm
    nwd = hd["nw"]["$h"]
    Yh  = HPM._harmonic_admittance(nwd, bus_idx)
    Yr  = Yh[non_ref, non_ref]
    Z   = Yr \ Matrix{ComplexF64}(I, length(non_ref), length(non_ref))
    for db in unit_buses
        nb = bus_idx[db]
        (nb in ref_ids) && continue
        z = Z[pos[nstar], pos[db]]
        argz = rad2deg(angle(z))
        vg   = seqsign(h) * (vg_ang[db] - vg_ang[8])   # relative VG rotation d→8
        imp  = argz - vg
        # wrap imp to (−180,180]
        imp = mod(imp + 180, 360) - 180
        @printf("%-4d %-8d %-10.4e %-10.2f %-10.2f %-10.2f\n", h, db, abs(z), argz, vg, imp)
    end
    println()
end
