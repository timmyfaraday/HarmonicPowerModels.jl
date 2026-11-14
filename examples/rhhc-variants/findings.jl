################################################################################
# examples/rhhc-variants/findings.jl                                           #
# Three findings of the RHHC paper on `ieee30-full-exact` (23 coordinates,     #
# harmonics 5–25, exact moments, Gaussian factor):                              #
#   1. the sign-predicted corners: the method without them, and its solution  #
#      checked at the corners and at independent random points;               #
#   2. the worst predicted corner of the binding constraint, in physical      #
#      inputs;                                                                 #
#   3. the √h scaling of the 50 Hz parameters against the physical model,     #
#      nominal and robust over factor intervals derived from the same box.    #
#                                                                               #
# SHHC_HMAX=25 julia -t 10 --project=examples/rhhc-variants \                   #
#     examples/rhhc-variants/findings.jl                                        #
# Output: output/ieee30-full-exact-findings.txt                                 #
################################################################################

isdefined(@__MODULE__, :STUDIES) || include(joinpath(@__DIR__, "builders.jl"))
using Printf, Random

const OUT = joinpath(STUDIES, "output", "ieee30-full-exact-findings.txt")
io = open(OUT, "w")
say(args...) = (println(io, args...); println(args...); flush(io))

s = study("cases/ieee30-full-exact.toml")
header(io, s)
d = length(s.box)
total(sol) = sum(values(sol))
violations(hd, g, ξ) = HPM._violations_at(hd, s.model, s.box, ξ, g)

# ── 1. with and without the sign-predicted corners ───────────────────────────
hd1 = deepcopy(s.hd)
t1 = @elapsed res = solve_rhhc(hd1, s.model, s.box, s.conic; pce_params = deepcopy(s.pce),
    hpf_optimizer = s.ipopt, posterior = 1024)
hd0 = deepcopy(s.hd)
t0 = @elapsed res0 = solve_rhhc(hd0, s.model, s.box, s.conic; pce_params = deepcopy(s.pce),
    hpf_optimizer = s.ipopt, posterior = 1024, corners = false)
dict    = res["dictionary"]
pos     = dict.limit .> 0.0
corners = dict.points[dict.labels .== :corner]
say("1. sign-predicted corners (", length(corners), " in the dictionary)")
for (name, r, t) in (("with corners", res, t1), ("without corners", res0, t0))
    g = Dict(h => I^2 for (h, I) in r["solution"])
    worst = maximum(maximum(violations(hd1, g, ξ)[pos]) for ξ in corners)
    say(@sprintf("   %-16s robust %.6f, %d iterations, %s, %.0f s; posterior worst %+.2e (%.1f %% violated); worst at the corners %+.2e",
        name, res === r ? r["lower_bound"] : r["lower_bound"], r["iterations"],
        r["signed_off"] ? "signed off" : "NOT signed off", t, r["posterior"].max_violation,
        100r["posterior"].violated_fraction, worst))
end

# ── 2. the worst predicted corner of the binding constraint ──────────────────
g  = Dict(h => I^2 for (h, I) in res["solution"])
V  = dictionary_violations(dict, g)
Vp = ifelse.(pos, V, -Inf)
j  = argmax([maximum(view(Vp, c, :)) for c in axes(Vp, 1)])   # binding constraint
ci = findall(==(:corner), dict.labels)
k  = ci[argmax(V[j, ci])]
p  = latent(s.box, dict.points[k])
say("2. binding constraint ", dict.keys[j], ": violation ", @sprintf("%+.2e", maximum(Vp[j, :])),
    " (", dict.labels[argmax(Vp[j, :])], "); its worst predicted corner ",
    @sprintf("%+.2e", V[j, k]), ":")
for (name, lo, hi) in zip(s.box.names, s.box.lower, s.box.upper)
    v = getfield(p, name)
    tag = isapprox(v, lo; rtol = 1e-9) ? "lower" : isapprox(v, hi; rtol = 1e-9) ? "upper" : "nominal"
    say(@sprintf("   %-22s %12.6g  (%s; interval %.6g to %.6g)", name, v, tag, lo, hi))
end

# ── 3. √h scaling of the 50 Hz parameters ────────────────────────────────────
cls = Dict(b => (s.details[b].kv < 100 ? :cable : :line) for b in keys(s.details))
sc  = ScaledBranchModel(s.hd, cls)                      # fundamental data = LCM, nominal
box_nom = ParameterBox(:kR_cable => (1.0, 1.0, 1.0))
hds = deepcopy(s.hd)
nom_s = solve_rhhc(hds, sc, box_nom, s.conic; pce_params = deepcopy(s.pce),
    hpf_optimizer = s.ipopt, nbackground = 0, use_seeds = false)
box0 = ParameterBox(s.box.names, s.box.nominal, s.box.nominal, s.box.nominal, s.box.scale)
nom_l = solve_rhhc(deepcopy(s.hd), s.model, box0, s.conic; pce_params = deepcopy(s.pce),
    hpf_optimizer = s.ipopt, nbackground = 0, use_seeds = false)
gs = Dict(h => I^2 for (h, I) in nom_s["master_solution"])
say("3. √h scaling of the 50 Hz parameters")
say(@sprintf("   nominal: scaled %.6f, physical %.6f (%+.1f %%); scaled solution on the physical model: violation %+.2e at the nominal point",
    nom_s["upper_bound"], nom_l["upper_bound"], 100(nom_s["upper_bound"] / nom_l["upper_bound"] - 1),
    maximum(violations(hd1, gs, zeros(d))[pos])))
# factor intervals per class from the 50 Hz data over the box: the corners and 1024 random points
Random.seed!(3)
pts = vcat(corners, [2 .* rand(d) .- 1 for _ in 1:1024])
ratio = Dict((c, q) => Float64[] for c in (:cable, :line), q in (:R, :L, :C))
for ξ in pts, (b, c) in cls
    x  = branch_data(s.model, b, latent(s.box, ξ), 1.0, :positive)
    x0 = sc.base[b]
    push!(ratio[(c, :R)], x.br_r / x0.br_r); push!(ratio[(c, :L)], x.br_x / x0.br_x)
    push!(ratio[(c, :C)], (x.b_fr + x.b_to) / (x0.b_fr + x0.b_to))
end
fbox = ParameterBox([Symbol("k", q, "_", c) => (1.0, minimum(ratio[(c, q)]), maximum(ratio[(c, q)]))
                     for c in (:cable, :line) for q in (:R, :L, :C)]...)
for (n, lo, hi) in zip(fbox.names, fbox.lower, fbox.upper)
    say(@sprintf("   factor %-9s %.4f to %.4f", n, lo, hi))
end
rob_s = solve_rhhc(deepcopy(s.hd), sc, fbox, s.conic; pce_params = deepcopy(s.pce),
    hpf_optimizer = s.ipopt, posterior = 256)
pc = posterior_check(hd1, s.model, s.box, rob_s["solution"]; n = 1024)
say(@sprintf("   robust over the factor box: %.6f (%s); on the physical model: worst %+.2e, violated at %.1f %% of 1024 random points (physical robust %.6f)",
    rob_s["lower_bound"], rob_s["signed_off"] ? "signed off" : "NOT signed off",
    pc.max_violation, 100pc.violated_fraction, res["lower_bound"]))
close(io)
