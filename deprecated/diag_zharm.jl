using HarmonicPowerModels, PowerModels, Ipopt, Printf
using LinearAlgebra
const PMs = PowerModels; const HPM = HarmonicPowerModels

include("examples/shhc/pce_params.jl")
data  = PMs.parse_file("examples/shhc/data/industrial_network_shhc.m")
H     = [1,5,7,11,13]
hdata = HPM.replicate(data, H=H)
HPM.update_hdata_with_pce_data!(hdata, PCE_PARAMS)
HPM.update_hdata_with_fundamental_hpf_results!(hdata, dHHC_NLP, Ipopt.Optimizer)

println("\n=== xfmr data at h=5 ===")
for (xid, xfmr) in sort(collect(hdata["nw"]["5"]["xfmr"]), by=x->parse(Int,x[1]))
    @printf("  xfmr %s: %d→%d  vg=%-6s  tr=%+.4f  ti=%+.4f  r1=%.5f  r2=%.5f  xsc=%.5f\n",
        xid, xfmr["f_bus"], xfmr["t_bus"], get(xfmr,"vg","?"),
        get(xfmr,"tr",0.0), get(xfmr,"ti",0.0),
        get(xfmr,"r1",0.0), get(xfmr,"r2",0.0), get(xfmr,"xsc",0.0))
end

fund_nw = hdata["nw"]["1"]
bus_ids_sorted = sort([bus["bus_i"] for (_,bus) in fund_nw["bus"]])
bus_idx = Dict(id => i for (i,id) in enumerate(bus_ids_sorted))

nw5 = hdata["nw"]["5"]
Yh = HPM._harmonic_admittance(nw5, bus_idx)

println("\nY_harm diagonal (bus_orig → matrix_idx):")
for (orig, idx) in sort(collect(bus_idx), by=x->x[1])
    @printf("  bus %d (idx %d): Y_diag = %+.4e %+.4ej\n", orig, idx, real(Yh[idx,idx]), imag(Yh[idx,idx]))
end

ref_bus_ids = Set(bus_idx[bus["bus_i"]] for (_,bus) in fund_nw["bus"] if bus["bus_type"]==3)
non_ref_ids  = [i for i in 1:length(bus_ids_sorted) if i ∉ ref_bus_ids]
non_ref_pos  = Dict(id => idx for (idx,id) in enumerate(non_ref_ids))
println("\nref_bus matrix indices: ", collect(ref_bus_ids))

Y_red = Yh[non_ref_ids, non_ref_ids]
nr    = length(non_ref_ids)
Z_red = Y_red \ Matrix{ComplexF64}(I, nr, nr)

n_pos = non_ref_pos[bus_idx[8]]
load_buses = [bus_idx[nw5["load"][d]["load_bus"]] for d in keys(nw5["load"])]
load_bus_positions = sort([non_ref_pos[b] for b in load_buses if b ∉ ref_bus_ids])
println("load bus positions in Z_red: ", load_bus_positions, " (for original buses 2,4,6,8)")

pce = hdata["pce"]
cp  = pce["cos_phi"][5][1]
sp  = pce["sin_phi"][5][1]
println("\nPCE trig at h=5: cos_phi=", round.(cp,digits=5), "  sin_phi=", round.(sp,digits=5))

println("\nIndividual Z_{bus8, d} (load d at position d_pos):")
Z_sum = 0.0 + 0.0im
for d_pos in load_bus_positions
    z = Z_red[n_pos, d_pos]
    @printf("  Z[bus8, pos%d] = %+.6e %+.6ej  |Z|=%.6e\n", d_pos, real(z), imag(z), abs(z))
    global Z_sum += z
end
@printf("Z_total = %+.6e %+.6ej  |Z|=%.6e\n", real(Z_sum), imag(Z_sum), abs(Z_sum))

alpha_k0 = Z_sum * (cp[1] + im*sp[1])
@printf("\nalpha[bus8,k=0] = %+.6e %+.6ej  |alpha|=%.6e\n", real(alpha_k0), imag(alpha_k0), abs(alpha_k0))
println("  (LP report says: ar=-23.07e+0  ai=+84.26e+0  from the txt file)")

# Check deterministic prediction
println("\n--- Check against deterministic dHHC_SOC ---")
# det: I_h=6.171e-4, injection angle=60° (crd=3.086e-4, cid=5.345e-4)
I_det = 6.1713e-4
ang_det = atan(5.3445e-4, 3.0857e-4)  # injection angle in radians
@printf("det I_h = %.4e pu,  angle = %.2f°\n", I_det, rad2deg(ang_det))
V8_det_from_Zharm = Z_sum * I_det * (cos(ang_det) + im*sin(ang_det))
@printf("LP prediction  |V8| = %.4e pu  (det Z model)\n", abs(V8_det_from_Zharm))
println("det actual     |V8| = 5.4601e-02 pu")
