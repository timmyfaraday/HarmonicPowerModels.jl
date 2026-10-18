using HarmonicPowerModels, PowerModels, Ipopt, Printf
using LinearAlgebra
const PMs = PowerModels; const HPM = HarmonicPowerModels

include("examples/shhc/pce_params.jl")
data  = PMs.parse_file("examples/shhc/data/industrial_network_shhc.m")
H     = [1,5,7,11,13]
hdata = HPM.replicate(data, H=H)
HPM.update_hdata_with_pce_data!(hdata, PCE_PARAMS)
HPM.update_hdata_with_fundamental_hpf_results!(hdata, dHHC_NLP, Ipopt.Optimizer)

# ── Build LP's Z_harm for h=5 and compute alpha per bus ──────────────────────
fund_nw = hdata["nw"]["1"]
bus_ids_sorted = sort([bus["bus_i"] for (_,bus) in fund_nw["bus"]])
bus_idx   = Dict(id => i for (i,id) in enumerate(bus_ids_sorted))
ref_ids   = Set(bus_idx[bus["bus_i"]] for (_,bus) in fund_nw["bus"] if bus["bus_type"]==3)
non_ref   = [i for i in 1:length(bus_ids_sorted) if i ∉ ref_ids]
non_ref_pos = Dict(id => idx for (idx,id) in enumerate(non_ref))

nw5   = hdata["nw"]["5"]
Yh    = HPM._harmonic_admittance(nw5, bus_idx)
Z_red = (Yh[non_ref, non_ref]) \ Matrix{ComplexF64}(I, length(non_ref), length(non_ref))

pce  = hdata["pce"]
cp0  = pce["cos_phi"][5][1][1]   # k=0
sp0  = pce["sin_phi"][5][1][1]

load_bus_positions = sort([non_ref_pos[bus_idx[nw5["load"][d]["load_bus"]]]
                           for d in keys(nw5["load"]) if bus_idx[nw5["load"][d]["load_bus"]] ∉ ref_ids])

# QCQP k=0 voltages at h=5 (from sHHC_QCQP results, cmd=1.36699e-3 per load)
cmd_qcqp = 1.36699e-3
qcqp_vr = Dict(0=>0.0, 1=>-2.8858e-4, 2=>2.3205e-3, 3=>2.3116e-3,
               4=>1.2659e-2, 5=>1.5297e-2, 6=>-2.7826e-3, 7=>-3.9047e-3, 8=>1.2284e-2)
qcqp_vi = Dict(0=>0.0, 1=>-1.6762e-3, 2=>-6.1380e-3, 3=>-6.5599e-3,
               4=>4.8767e-2, 5=>5.1542e-2, 6=>2.1365e-2, 7=>1.2696e-2, 8=>-5.1510e-2)

println("\nBus-by-bus comparison  (LP alpha[n,k=0] * cmd_qcqp  vs  QCQP v[n,k=0])")
println("="^80)
@printf("%-6s  %-22s  %-22s  %-10s  %-10s\n",
        "Bus", "LP prediction Vr+jVi", "QCQP actual Vr+jVi", "|LP|/cmd", "|QCQP|/cmd")
println("-"^80)
for bus_orig in sort(collect(keys(qcqp_vr)))
    if bus_orig == 0
        @printf("%-6d  %+.4e%+.4ej  %+.4e%+.4ej  %8s  %8s\n",
                bus_orig, 0.0, 0.0, 0.0, 0.0, "(ref)", "(ref)")
        continue
    end
    n_pos = non_ref_pos[bus_idx[bus_orig]]
    Zt    = sum(Z_red[n_pos, d_pos] for d_pos in load_bus_positions)
    alpha = Zt * (cp0 + im*sp0)
    V_LP  = alpha * cmd_qcqp
    V_QCQP = qcqp_vr[bus_orig] + im*qcqp_vi[bus_orig]
    @printf("%-6d  %+.4e%+.4ej  %+.4e%+.4ej  %10.2f  %10.2f\n",
            bus_orig, real(V_LP), imag(V_LP), real(V_QCQP), imag(V_QCQP),
            abs(Zt)*abs(cp0+im*sp0), abs(V_QCQP)/cmd_qcqp)
end

# ── Verify individual transformer equations with QCQP k=0 voltages ───────────
println("\n=== Transformer equation check at h=5, k=0 ===")
for (xid, xfmr) in sort(collect(nw5["xfmr"]), by=x->parse(Int,x[1]))
    f = xfmr["f_bus"]; t = xfmr["t_bus"]
    tr = xfmr["tr"];   ti = xfmr["ti"]
    r1 = xfmr["r1"];   r2 = xfmr["r2"];  xsc = xfmr["xsc"]
    vg = get(xfmr,"vg","?")

    Vf = qcqp_vr[f] + im*qcqp_vi[f]
    Vt = qcqp_vr[t] + im*qcqp_vi[t]
    a  = tr + im*ti

    # Series model: I_f = (V_f - a*V_t) / (r1+r2+j*xsc)
    z_series = (r1+r2) + im*xsc
    I_f_series = (Vf - a*Vt) / z_series
    I_t_series = -conj(a)*I_f_series

    # How much does the LP expect as current through this xfmr?
    # From KCL at t-bus: I entering t = sum of load injections downstream
    # Downstream loads of each xfmr:
    downstream_loads = Dict(
        "1"=>[1,2,3,4], "2"=>[2,3,4], "3"=>[3,4], "4"=>[4]
    )
    dl = get(downstream_loads, xid, Int[])
    I_expected = length(dl) * cmd_qcqp * (cp0 + im*sp0)

    @printf("\nxfmr %s (%s %d→%d)  a=%+.4f%+.4fj  z_series=%+.4f%+.4fj\n",
            xid, vg, f, t, tr, ti, r1+r2, xsc)
    @printf("  V_f = %+.4e%+.4ej   V_t = %+.4e%+.4ej\n", real(Vf), imag(Vf), real(Vt), imag(Vt))
    @printf("  I_f (LP series model) = %+.4e%+.4ej  |I_f|=%.4e\n", real(I_f_series), imag(I_f_series), abs(I_f_series))
    @printf("  I_expected (from %d downstream loads) = %+.4e%+.4ej  |I|=%.4e\n",
            length(dl), real(I_expected), imag(I_expected), abs(I_expected))
    ratio = abs(I_f_series) > 1e-15 ? abs(I_expected)/abs(I_f_series) : Inf
    @printf("  |I_expected| / |I_f_LP| = %.4f\n", ratio)
end
