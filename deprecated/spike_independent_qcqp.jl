# spike_independent_qcqp.jl  —  R1 feasibility spike (THROWAWAY prototype)
#
# Question: does the existing K-parametric sHHC_QCQP code solve when fed
# MULTIVARIATE independent-germ PCE data, and how does Ipopt scale with
# P = C(N+δ, δ)?  No production refactor — we replicate the solve_shhc(sHHC_QCQP)
# pipeline but inject a multivariate independent hdata["pce"] (one Beta germ per
# harmonic unit) instead of the shared single-germ data.
#
# Run:  julia --project spike_independent_qcqp.jl

using HarmonicPowerModels, PowerModels
using Ipopt
import JuMP
using Printf

const PMs = PowerModels
const HPM = HarmonicPowerModels
include(joinpath(HPM.BASE_DIR, "examples", "shhc", "pce_params.jl"))
const SILENT = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes")

# ── multivariate independent data builder (the prototype of Phase A) ──────────
# pce_params may carry "harmonic_units" => Vector{Int} (load ids that are
# harmonic sources).  Default = all loads.  Normal loads get zero cos/sin.
function multivariate_pce!(hdata, pce_params)
    α̂  = pce_params["alpha_hat"]; β̂ = pce_params["beta_hat"]
    δ  = pce_params["delta"];     nt = get(pce_params, "n_terms", 4)
    ε  = pce_params["epsilon"];   hb = pce_params["harmonics"]

    all_loads = sort(unique(parse(Int, d)
                            for (nw, ntw) in hdata["nw"] if nw != "1"
                            for d in keys(ntw["load"])))
    hunits = sort(get(pce_params, "harmonic_units", all_loads))
    N      = length(hunits)
    dimmap = Dict(uid => i for (i, uid) in enumerate(hunits))

    pb = HPM.PCEBasis(N, α̂, β̂, δ)
    P  = pb.dim
    M  = HPM.dense_tensor(pb)

    cos_phi = Dict{Int, Dict{Int, Vector{Float64}}}()
    sin_phi = Dict{Int, Dict{Int, Vector{Float64}}}()
    for (nw_str, ntw) in hdata["nw"]
        nw = parse(Int, nw_str)
        (nw == 1 || !haskey(hb, nw)) && continue
        vmin, vmax = hb[nw]
        cos_phi[nw] = Dict{Int, Vector{Float64}}()
        sin_phi[nw] = Dict{Int, Vector{Float64}}()
        for d_str in keys(ntw["load"])
            d = parse(Int, d_str)
            if haskey(dimmap, d)
                cc, sc = HPM.rotated_trig_pce(pb, dimmap[d], vmin, vmax; n_terms = nt)
            else
                cc, sc = zeros(P), zeros(P)
            end
            cos_phi[nw][d] = cc
            sin_phi[nw][d] = sc
        end
    end

    hdata["pce"] = Dict("basis" => pb, "M" => M, "lambda" => HPM.cantelli_lambda(ε),
                        "psi_norms" => pb.sp2, "K" => 0:(P-1),
                        "cos_phi" => cos_phi, "sin_phi" => sin_phi)
    return hdata
end

# replicate the body of solve_shhc(sHHC_QCQP) but with an explicit pce builder
function solve_qcqp_with(hdata, params, builder; optimizer=SILENT)
    builder(hdata, params)
    HPM.update_hdata_with_fundamental_hpf_results!(hdata, dHHC_NLP, optimizer)
    HPM.update_hdata_with_fairness_principle_data!(hdata, sHHC_QCQP, optimizer)
    return PMs.solve_model(hdata, sHHC_QCQP, optimizer, HPM.build_shhc;
        ref_extensions      = [HPM.ref_add_filter!, HPM.ref_add_xfmr!],
        solution_processors = [HPM.sol_data_model!],
        multinetwork        = true)
end

get_I(sol, h) = (s = sol["solution"]; haskey(s, "nw") ? s["nw"]["$h"]["load"]["1"]["cmd"] : s[h])
status_of(sol) = string(sol["termination_status"])

# ══════════════════════════════════════════════════════════════════════════════
# Part 1 — industrial network, N=4, correctness + directional check
# ══════════════════════════════════════════════════════════════════════════════
println("="^70)
println("PART 1 — industrial network (N=4 harmonic units, P=15)")
println("="^70)

ind_path = joinpath(HPM.BASE_DIR, "examples", "shhc", "data", "industrial_network_shhc.m")
ind_data = PMs.parse_file(ind_path)
H    = [1, 5, 7, 11, 13]
harm = [5, 7, 11, 13]

t_sh = @elapsed sol_sh = solve_qcqp_with(HPM.replicate(ind_data, H=H), PCE_PARAMS,
                                         HPM.update_hdata_with_pce_data!)
t_in = @elapsed sol_in = solve_qcqp_with(HPM.replicate(ind_data, H=H), PCE_PARAMS,
                                         multivariate_pce!)

# independent LP reduction (IHD-only ground reference)
beta_params = (alpha_shape=PCE_PARAMS["alpha_hat"], beta_shape=PCE_PARAMS["beta_hat"], delta=PCE_PARAMS["delta"])
lp_ind = Dict{Int,Any}()
for h in harm
    hbv  = PCE_PARAMS["harmonics"][h]
    uids = sort(parse.(Int, collect(keys(ind_data["load"]))))
    ap   = Dict(u => (rad2deg(hbv.v_min), rad2deg(hbv.v_max)) for u in uids)
    lp_ind[h] = HPM.compare_shared_vs_independent(ind_data, h, PCE_PARAMS["epsilon"], beta_params, ap)
end

println("\n  shared QCQP      : $(status_of(sol_sh))  ($(round(t_sh, digits=2)) s)")
println("  independent QCQP : $(status_of(sol_in))  ($(round(t_in, digits=2)) s)")
@printf("\n  %-4s %-15s %-15s %-15s %-15s\n", "h", "I_QCQP_shared", "I_QCQP_indep", "I_LP_indep", "indep<shared?")
for h in harm
    Ish = get_I(sol_sh, h); Iin = get_I(sol_in, h); Ilp = lp_ind[h].I_star_independent
    @printf("  %-4d %-15.6e %-15.6e %-15.6e %-15s\n", h, Ish, Iin, Ilp, Iin < Ish*0.999 ? "yes" : "no/≈")
end
println("\n  (independent germs de-correlate angles → more variance → I_indep < I_shared expected)")

# ══════════════════════════════════════════════════════════════════════════════
# Part 2 — IEEE-30 scaling sweep: solve time vs N (hence P)
# ══════════════════════════════════════════════════════════════════════════════
println("\n" * "="^70)
println("PART 2 — IEEE-30 scaling: solve time vs N harmonic units")
println("="^70)

function run_scaling()
    ieee_path = joinpath(HPM.BASE_DIR, "examples", "shhc", "data", "30bus_network_shhc.m")
    ieee_data = PMs.parse_file(ieee_path)
    all_uids  = sort(parse.(Int, collect(keys(ieee_data["load"]))))
    println("  IEEE-30 has $(length(all_uids)) loads total\n")

    # baseline: does the SHARED QCQP even solve IEEE-30 in this setup?
    t_base, base_status = NaN, "?"
    try
        t_base = @elapsed sol_base = solve_qcqp_with(HPM.replicate(ieee_data, H=H), PCE_PARAMS,
                                                     HPM.update_hdata_with_pce_data!)
        base_status = status_of(sol_base)
    catch e
        base_status = "ERROR: " * first(split(sprint(showerror, e), '\n'))
    end
    @printf("  shared baseline (N=1, P=3): %-22s (%s s)\n\n", base_status, isnan(t_base) ? "—" : round(t_base, digits=2))

    @printf("  %-4s %-6s %-12s %-22s\n", "N", "P", "solve [s]", "status")
    for N in [2, 4, 6, 8, 11, 15, 21]
        N > length(all_uids) && continue
        P = binomial(N + PCE_PARAMS["delta"], PCE_PARAMS["delta"])
        params = merge(PCE_PARAMS, Dict("harmonic_units" => all_uids[1:N]))
        st, status = NaN, "?"
        try
            st = @elapsed sol = solve_qcqp_with(HPM.replicate(ieee_data, H=H), params, multivariate_pce!)
            status = status_of(sol)
        catch e
            status = "ERROR: " * first(split(sprint(showerror, e), '\n'))
        end
        @printf("  %-4d %-6d %-12s %-22s\n", N, P, isnan(st) ? "—" : string(round(st, digits=2)), status)
        flush(stdout)
    end
end

run_scaling()
println("\nSpike complete.")
