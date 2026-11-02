################################################################################
# examples/rhhc/ieee30_robust.jl                                               #
# Robust harmonic hosting capacity on the IEEE 30-bus system with 33 kV       #
# cables (scenario `inductive_cable33` of the SHHC paper), 21 harmonic units, #
# local angle reference, ε = 0.05.                                             #
#                                                                               #
#   1. nominal stochastic (sHHC_SOC) versus robust (solve_rhhc), Cantelli and  #
#      Gaussian bounds, each checked a posteriori on 1024 random points;       #
#   2. ablation of the oracle's dictionary: with and without the modal seeds,  #
#      for a low-dimensional box (one factor per quantity and branch class,    #
#      d = 5, vertices included) and a high-dimensional one (a length factor   #
#      per cable, d = 25, no vertices), against the size of the background.    #
#                                                                               #
# The network model is analytic (ScaledBranchModel); the cable and line data  #
# of LineCableModels.jl replace it once the IEEE30 mapping is available.      #
# Parts: `julia ieee30_robust.jl [nominal] [d5] [d25]` (default: all three).  #
# Output: examples/rhhc/output/ieee30_robust_<parts>.txt                      #
################################################################################

using HarmonicPowerModels, PowerModels, JuMP, Ipopt, MosekTools, Printf
const HPM = HarmonicPowerModels
PowerModels.silence()

const EX = joinpath(@__DIR__, "..", "shhc")
include(joinpath(EX, "pce_params.jl"))
include(joinpath(EX, "ieee30_scenarios.jl"))

ipopt = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
mosek = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)
harm  = [5, 7, 11, 13]; H = [1; harm]
rng   = Dict(h => (PCE_PARAMS["harmonics"][h].v_min, PCE_PARAMS["harmonics"][h].v_max) for h in harm)

"Replicated data with `units` as the only harmonic units, local angle reference."
function restricted(data, units)
    hd = HPM.replicate(deepcopy(data), H = H); hd["angle_reference"] = "local"
    for h in harm, k in collect(keys(hd["nw"]["$h"]["load"]))
        parse(Int, k) in units || delete!(hd["nw"]["$h"]["load"], k)
    end
    return hd
end
function pce(units, bound)
    p = Dict{String,Any}(deepcopy(PCE_PARAMS)...)
    p["harmonics"] = Dict{Int,Any}(h => rng[h] for h in harm)
    p["harmonic_units"] = units; p["bound"] = bound
    return p
end

data  = ieee30_scenario("inductive_cable33", ipopt)
units = sort(parse.(Int, collect(keys(data["load"]))))[1:21]
isxfmr(b) = get(b, "transformer", false) || !(get(b, "tap", 1.0) in (0.0, 1.0)) ||
            get(b, "shift", 0.0) != 0.0
cables = sort([parse(Int, k) for (k, b) in data["branch"] if get(b, "b_fr", 0.0) > 0 && !isxfmr(b)])
ohls   = sort([parse(Int, k) for (k, b) in data["branch"] if get(b, "b_fr", 0.0) == 0 && !isxfmr(b)])
classes = merge(Dict(b => :cable for b in cables), Dict(b => :ohl for b in ohls))

box5 = ParameterBox(:kC_cable => (1.0, 0.90, 1.10), :kL_cable => (1.0, 0.95, 1.05),
                    :kR_cable => (1.0, 0.80, 1.20), :kL_ohl   => (1.0, 0.95, 1.05),
                    :kR_ohl   => (1.0, 0.80, 1.20))
box25 = ParameterBox(:kC_cable => (1.0, 0.90, 1.10), :kL_cable => (1.0, 0.95, 1.05),
                     :kR_cable => (1.0, 0.80, 1.20),
                     [Symbol("len_", b) => (1.0, 0.95, 1.05) for b in cables]...)

parts  = isempty(ARGS) ? ["nominal", "d5", "d25"] : ARGS
outdir = joinpath(@__DIR__, "output"); mkpath(outdir)
io = open(joinpath(outdir, "ieee30_robust_" * join(parts, "_") * ".txt"), "w")
out(args...) = (println(io, args...); println(args...); flush(io))
out("IEEE 30-bus, inductive_cable33, 21 units, ε = 0.05; ", length(cables), " cables, ",
    length(ohls), " overhead lines (one class each)")

# ── 1. nominal versus robust ──────────────────────────────────────────────────
"nominal" in parts && out("\n== 1. nominal stochastic versus robust, box d = 5 ==")
for bound in ("cantelli", "gaussian")
    "nominal" in parts || break
    hd  = restricted(data, units)
    nom = HPM.solve_shhc(hd, sHHC_SOC, mosek; pce_params = pce(units, bound), hpf_optimizer = ipopt)
    model = ScaledBranchModel(hd, classes)
    pn  = posterior_check(hd, model, box5, nom["solution"]; n = 1024)
    hd  = restricted(data, units)
    model = ScaledBranchModel(hd, classes)
    t = @elapsed res = solve_rhhc(hd, model, box5, mosek; pce_params = pce(units, bound),
                                  hpf_optimizer = ipopt, posterior = 1024, verbose = true)
    pr = res["posterior"]
    out(@sprintf("%-8s nominal %.6f (violated at %.1f %% of the box, worst %+.3f)", bound,
                 nom["objective_value"], 100pn.violated_fraction, pn.max_violation))
    out(@sprintf("%-8s robust  %.6f  [LB %.6f, UB %.6f, gap %.1e]  %d iterations, %.0f s, %s",
                 bound, res["lower_bound"], res["lower_bound"], res["upper_bound"], res["gap"],
                 res["iterations"], t, res["signed_off"] ? "signed off" : "NOT signed off"))
    out(@sprintf("%-8s robust  a posteriori: worst %+.2e (%s), violated at %.1f %%", bound,
                 pr.max_violation, string(pr.constraint), 100pr.violated_fraction))
    out("         per harmonic, robust / nominal: ",
        join([@sprintf("h%d %.3f", h, res["solution"][h] / nom["solution"][h]) for h in harm], ", "))
    out("         history (v, UB, cuts): ",
        join([@sprintf("(%.3g, %.5f, %d)", r.max_violation, r.upper_bound, r.cuts) for r in res["history"]], " "))
end

# ── 2. ablation ───────────────────────────────────────────────────────────────
ablations = [(n, b, bg) for (n, b, bg) in (("d5", box5, (0, 64, 1024)), ("d25", box25, (64, 512)))
             if n in parts]
isempty(ablations) || out("\n== 2. ablation of the dictionary (Cantelli) ==")
isempty(ablations) || out("box     seeds  background  iterations  robust(LB)   gap      time[s]  posterior worst  violated")
for (name, box, bgs) in ablations
    for nb in bgs, use_seeds in (false, true)
        hd    = restricted(data, units)
        model = ScaledBranchModel(hd, classes; lengths = cables)
        t = @elapsed res = solve_rhhc(hd, model, box, mosek; pce_params = pce(units, "cantelli"),
                                      hpf_optimizer = ipopt, nbackground = nb, use_seeds,
                                      posterior = 1024, verbose = true)
        pr = res["posterior"]
        out(@sprintf("%-6s  %-5s  %10d  %10d  %.6f  %.1e  %7.0f  %+15.2e  %7.1f %%", name,
                     use_seeds ? "yes" : "no", nb, res["iterations"], res["lower_bound"],
                     res["gap"], t, pr.max_violation, 100pr.violated_fraction))
    end
end
close(io)
