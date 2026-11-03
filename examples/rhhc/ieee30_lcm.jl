################################################################################
# examples/rhhc/ieee30_lcm.jl                                                  #
# Starter: LineCableModels.jl cables in the IEEE 30-bus robust test case.      #
#                                                                               #
# The 22 lines of the 33 kV level are underground cables (scenario             #
# `inductive_cable33` of examples/shhc). Here their per-unit-length data come  #
# from LineCableModels: one cable design shared by all 22, built from a named  #
# parameter tuple p, and the exact π-model of each cable at its own length.    #
# The 132 kV overhead lines and the transformers keep the case data.           #
#                                                                               #
# To plug in your own work, edit the two marked blocks:                        #
#   (1) `cable_zy`  — the cable design and system; return the positive-        #
#                     sequence z1 [Ω/m] and y1 [S/m] at frequency f for p;     #
#   (2) `box`       — the uncertain inputs and their bounds; every name in the #
#                     box must be a field that `cable_zy` reads from p.        #
#                                                                               #
# Run:  julia -t 10 --project=<env> ieee30_lcm.jl                              #
# Env:  HarmonicPowerModels (branch feature/rhhc-cutting-plane, developed),    #
#       LineCableModels (main, developed), PowerModels, JuMP, Ipopt, and a     #
#       conic solver (MosekTools if licensed, else Clarabel).                  #
# Output: examples/rhhc/output/ieee30_lcm.txt                                  #
################################################################################

using HarmonicPowerModels, PowerModels, LineCableModels, JuMP, Ipopt, Printf, LinearAlgebra
const HPM = HarmonicPowerModels
PowerModels.silence()
BLAS.set_num_threads(1)     # small matrices: BLAS threads only compete with the oracle's tasks

const EX = joinpath(@__DIR__, "..", "shhc")
include(joinpath(EX, "pce_params.jl"))
include(joinpath(EX, "ieee30_scenarios.jl"))

ipopt = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
conic = try
    @eval using MosekTools
    optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)
catch
    @eval using Clarabel
    optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false)
end

# ── (1) the cable: positive-sequence z1, y1 per meter at frequency f ───────────
# One 18/30 kV XLPE cable type, 1x630 Al, trefoil at 1 m depth, sheaths grounded at both
# ends (Kron-eliminated), phases ideally transposed so z1 = zs − zm. Replace the design,
# the formation, the bonding and the earth with yours; keep the signature.
const CABLE_NOMINAL = (rho_core = 2.8264e-8, r_core = 14.16e-3, t_xlpe = 8.0e-3,
                       eps_xlpe = 2.3, rho_earth = 100.0)

function cable_zy(p, f)
    core_metal = Material(:conductor, p.rho_core, 1.0, 1.0, 20.0, 0.00403)
    semicon    = Material(:semicon, 1_000.0, 1_000.0)
    xlpe       = Material(:insulator, 1.0e14, p.eps_xlpe)
    copper     = Material(:conductor, 1.7241e-8, 1.0, 1.0, 20.0, 0.00393)
    pe         = Material(:insulator, 1.0e14, 2.3)
    design = @cable "xlpe_18_30kv_1x630al" begin
        @terminal :core begin
            core(core_metal; r = p.r_core, tag = :core)
        end
        screen(semicon; t = 0.8e-3, tag = :semicon_inner)
        insulation(xlpe; t = p.t_xlpe, tag = :insulation)
        screen(semicon; t = 0.8e-3, tag = :semicon_outer)
        @terminal :sheath begin
            sheath(copper; t = 0.3e-3, tag = :screen)
        end
        jacket(pe; t = 2.5e-3, tag = :jacket)
    end
    d = 2 * r_ex(design.outer)
    placements = [Placement(-d / 2, -1.0), Placement(d / 2, -1.0),
                  Placement(0.0, -1.0 + d * sqrt(3) / 2)]
    system = build(CableSystem, fill(design, 3), placements;
        connections = [(core = i, sheath = 0) for i in 1:3], id = "cable", line_length = 1.0)
    problem = LineParametersProblem(system, homogeneous(rho = p.rho_earth, eps_r = 10.0);
        temperature = 20.0, frequencies = [f])
    r = compute(problem, CoaxialFormulation(transposition = true))
    return r.Z[1, 1, 1] - r.Z[1, 2, 1], r.Y[1, 1, 1] - r.Y[1, 2, 1]
end

# ── (2) the parameter box: uncertain inputs and bounds ──────────────────────────
# name => (nominal, lower, upper[, :log]); a name not in the box keeps CABLE_NOMINAL.
box = ParameterBox(:eps_xlpe  => (2.3, 2.25, 2.5),
                   :t_xlpe    => (8.0e-3, 0.95 * 8.0e-3, 1.10 * 8.0e-3),
                   :r_core    => (14.16e-3, 0.99 * 14.16e-3, 1.01 * 14.16e-3),
                   :rho_earth => (100.0, 10.0, 1000.0, :log))

# ── nothing below needs editing ────────────────────────────────────────────────

# The oracle evaluates the network at many points p and frequencies f, from several threads:
# all 22 cables share one design, so one LineCableModels solve per (p, f) serves all of them.
const ZY_CACHE = Dict{Tuple{NamedTuple,Float64},Tuple{ComplexF64,ComplexF64}}()
const ZY_LOCK  = ReentrantLock()
function cable(p, f, seq)
    seq === :zero && error("no zero-sequence data yet (z0 = zs + 2zm, y0 = ys + 2ym)")
    q   = merge(CABLE_NOMINAL, p)
    key = (q, Float64(f))
    v = lock(() -> get(ZY_CACHE, key, nothing), ZY_LOCK)
    v === nothing || return v
    v = cable_zy(q, f)
    lock(() -> (ZY_CACHE[key] = v), ZY_LOCK)
    return v
end

harm = [5, 7, 11, 13]; H = [1; harm]
rng  = Dict(h => (PCE_PARAMS["harmonics"][h].v_min, PCE_PARAMS["harmonics"][h].v_max) for h in harm)

data  = ieee30_scenario("inductive_cable33", ipopt)
units = sort(parse.(Int, collect(keys(data["load"]))))[1:21]          # the 21 harmonic units
isxfmr(b) = get(b, "transformer", false) || !(get(b, "tap", 1.0) in (0.0, 1.0)) ||
            get(b, "shift", 0.0) != 0.0
cables = sort([parse(Int, k) for (k, b) in data["branch"] if get(b, "b_fr", 0.0) > 0 && !isxfmr(b)])

# each cable at its own route length (from the case), 33 kV on 100 MVA
zbase = 33.0^2 / data["baseMVA"]
model = ExactPiModel(Dict(b => (zy = cable, len = 1e3 * data["branch"]["$b"]["cable_km"], zbase = zbase)
                          for b in cables))

# the fundamental power flow uses the nominal LineCableModels cables as well
for b in cables, (k, v) in pairs(branch_data(model, b, nominal_point(box), 1.0, :positive))
    data["branch"]["$b"][string(k)] = v
end

function prepare(data)
    hd = HPM.replicate(deepcopy(data), H = H); hd["angle_reference"] = "local"
    for h in harm, k in collect(keys(hd["nw"]["$h"]["load"]))
        parse(Int, k) in units || delete!(hd["nw"]["$h"]["load"], k)
    end
    return hd
end
pce = Dict{String,Any}(deepcopy(PCE_PARAMS)...)
pce["harmonics"] = Dict{Int,Any}(h => rng[h] for h in harm)
pce["harmonic_units"] = units; pce["bound"] = "cantelli"

outdir = joinpath(@__DIR__, "output"); mkpath(outdir)
io = open(joinpath(outdir, "ieee30_lcm.txt"), "w")
out(args...) = (println(io, args...); println(args...); flush(io))
out("IEEE 30-bus, inductive_cable33, ", length(cables), " cables from LineCableModels, ",
    length(units), " units, ε = 0.05, Cantelli, ", Threads.nthreads(), " threads")

# nominal stochastic capacity with the LineCableModels cables: a zero-width box, so the
# cables are evaluated by LineCableModels at every harmonic frequency (solve_shhc on the
# replicated data would instead scale their 50 Hz data by r·√h, x·h, b·h)
box0 = ParameterBox(box.names, box.nominal, box.nominal, box.nominal, box.scale)
hd   = prepare(data)
nom  = solve_rhhc(hd, model, box0, conic; pce_params = deepcopy(pce), hpf_optimizer = ipopt,
                  nbackground = 0, use_seeds = false)
out(@sprintf("nominal stochastic: %.6f (%s)", nom["upper_bound"], nom["termination_status"]))
out("resonances at the nominal point: ",
    join([@sprintf("h = %.2f", r.h) for r in modal_scan(hd, model, nominal_point(box))], ", "))

# robust capacity over the box
hd = prepare(data)
t  = @elapsed res = solve_rhhc(hd, model, box, conic; pce_params = deepcopy(pce),
                               hpf_optimizer = ipopt, posterior = 1024, verbose = true)
pr = res["posterior"]
out(@sprintf("robust: %.6f  [LB %.6f, UB %.6f, gap %.1e], %d iterations, %.0f s, %s",
             res["lower_bound"], res["lower_bound"], res["upper_bound"], res["gap"],
             res["iterations"], t, res["signed_off"] ? "signed off" : "NOT signed off"))
out("per harmonic, robust / nominal: ",
    join([@sprintf("h%d %.3f", h, res["solution"][h] / nom["master_solution"][h]) for h in harm], ", "))
out(@sprintf("a posteriori (1024 random points): worst %+.2e (%s), violated at %.1f %%",
             pr.max_violation, string(pr.constraint), 100pr.violated_fraction))
out("LineCableModels solves: ", length(ZY_CACHE))
close(io)
