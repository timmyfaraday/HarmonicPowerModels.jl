################################################################################
# HarmonicPowerModels.jl — SHHC paper §IV-B: equivalence of the methods          #
#                                                                               #
# On the purely inductive IEEE 30-bus base network:                              #
#   1. point mass (all angle ranges collapsed to their mean, N = 4): dHHC_SOC,   #
#      sHHC_SOC and sHHC_QCQP give the same hosting capacity;                    #
#   2. stochastic (N = 4): sHHC_SOC and sHHC_QCQP agree for the Cantelli,        #
#      Gaussian and empirical risk factor;                                       #
#   3. reduction (N = 21): the hosting capacity reduced to the binding           #
#      individual harmonic voltage limits equals dHHC_SOC at the mean angles,   #
#      and its exact 95 % value equals sHHC_SOC with the empirical risk factor.  #
#                                                                               #
# Output: output/paper/equivalence.csv (check, method A, method B, A, B, Δ%)    #
# Author: Tom Van Acker                                                          #
################################################################################

include(joinpath(@__DIR__, "paper_common.jl"))

io = open(joinpath(PAPER_OUT, "equivalence.csv"), "w")
csvline(io, "check", "N", "method_a", "method_b", "value_a", "value_b", "rel_diff_pct")
row(check, N, a, b, va, vb) = csvline(io, check, N, a, b, @sprintf("%.8e", va), @sprintf("%.8e", vb),
                                      @sprintf("%.2e", 100 * (vb / va - 1)))

data = ieee30_scenario("inductive", ipopt)
Eφ(h) = rng[h][1] + (rng[h][2] - rng[h][1]) * Eω
dhhc_total(sol, u1) = sum(sol["solution"]["nw"]["$h"]["load"]["$u1"]["cmd"] for h in harm)

# ── 1. point mass, N = 4 ───────────────────────────────────────────────────────
u4 = paper_units(data, 4)
hd = restricted(data, u4)
for h in harm, d in u4; hd["nw"]["$h"]["load"]["$d"]["ref_angle"] = rad2deg(Eφ(h)); end
det = dhhc_total(HPM.solve_hhc(hd, dHHC_SOC, mosek, ipopt), u4[1])
pm  = paper_pce(u4, nothing); pm["harmonics"] = Dict{Int,Any}(h => (Eφ(h), Eφ(h)) for h in harm)
soc = sum(values(HPM.solve_shhc(restricted(data, u4), sHHC_SOC, mosek; pce_params = deepcopy(pm), hpf_optimizer = ipopt)["solution"]))
qc  = dhhc_total(HPM.solve_shhc(restricted(data, u4), sHHC_QCQP, ipopt; pce_params = deepcopy(pm)), u4[1])
row("point mass", 4, "dHHC_SOC", "sHHC_SOC", det, soc)
row("point mass", 4, "dHHC_SOC", "sHHC_QCQP", det, qc)

# ── 2. stochastic, N = 4 ───────────────────────────────────────────────────────
for b in ("cantelli", "gaussian")
    p = paper_pce(u4, nothing; bound = b)
    s = sum(values(HPM.solve_shhc(restricted(data, u4), sHHC_SOC, mosek; pce_params = deepcopy(p), hpf_optimizer = ipopt)["solution"]))
    q = dhhc_total(HPM.solve_shhc(restricted(data, u4), sHHC_QCQP, ipopt; pce_params = deepcopy(p)), u4[1])
    row("stochastic $b", 4, "sHHC_SOC", "sHHC_QCQP", s, q)
end
pe = paper_pce(u4, nothing)
se = sum(values(HPM.solve_shhc_calibrated(restricted(data, u4), sHHC_SOC, mosek; pce_params = deepcopy(pe), hpf_optimizer = ipopt, nsamp = 20000)["solution"]))
qe = HPM.solve_shhc_calibrated(restricted(data, u4), sHHC_QCQP, ipopt; pce_params = deepcopy(pe), nsamp = 20000)
row("stochastic empirical", 4, "sHHC_SOC", "sHHC_QCQP", se, dhhc_total(qe, u4[1]))

# ── 3. reduction, N = 21 ───────────────────────────────────────────────────────
u21 = paper_units(data)
r   = recombination("inductive", NS, "local"; data = data)
hd  = restricted(data, u21)
for h in harm, d in u21; hd["nw"]["$h"]["load"]["$d"]["ref_angle"] = rad2deg(Eφ(h)); end
row("reduction", NS, "dHHC_SOC", "reduced deterministic", dhhc_total(HPM.solve_hhc(hd, dHHC_SOC, mosek, ipopt), u21[1]), tot(r, :det))
row("reduction", NS, "sHHC_SOC empirical", "exact 95 % quantile", sum(values(soc_hosting(data, u21, nothing, "empirical"))), tot(r, :q95))
close(io)
println("✓ output/paper/equivalence.csv")
