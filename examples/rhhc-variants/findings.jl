# Evaluate three RHHC comparisons on ieee30-full-exact: 23 coordinates,
# characteristic harmonics from 5 to 25, exact moments and the Gaussian factor.
# Compare the cutting-set method with and without the predicted corners, checking
# both solutions at those corners and at independent random points. Report the
# physical inputs at the worst predicted corner of the binding constraint.
# Compare the physical model with the square-root harmonic scaling of its 50 Hz
# parameters: nominal capacity and robust hosting capacity over factor intervals
# derived from the physical box.
#
# Run with SHHC_HMAX=25 as shown in the README.
# The report is output/ieee30-full-exact-findings.txt.

isdefined(@__MODULE__, :STUDIES) || include(joinpath(@__DIR__, "builders.jl"))
using Printf

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
dict = res["dictionary"]
pos = dict.limit .> 0.0
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
g = Dict(h => I^2 for (h, I) in res["solution"])
V = dictionary_violations(dict, g)
Vp = ifelse.(pos, V, -Inf)
j = argmax([maximum(view(Vp, c, :)) for c in axes(Vp, 1)])   # binding constraint
ci = findall(==(:corner), dict.labels)
k = ci[argmax(V[j, ci])]
p = latent(s.box, dict.points[k])
say("2. binding constraint ", dict.keys[j],
    ": violation ", @sprintf("%+.2e", maximum(Vp[j, :])),
    " (", dict.labels[argmax(Vp[j, :])], "); its worst predicted corner ",
    @sprintf("%+.2e", V[j, k]), ":")
for (name, lo, hi) in zip(s.box.names, s.box.lower, s.box.upper)
    v = getfield(p, name)
    tag = isapprox(v, lo; rtol = 1e-9) ? "lower" :
          isapprox(v, hi; rtol = 1e-9) ? "upper" : "nominal"
    say(@sprintf("   %-22s %12.6g  (%s; interval %.6g to %.6g)", name, v, tag, lo, hi))
end

# ── 3. √h scaling of the 50 Hz parameters ────────────────────────────────────
scaled = study("cases/ieee30-full-exact-scaled.toml")
sc, fbox = scaled.model, scaled.box
box_nom = ParameterBox(fbox.names, fbox.nominal, fbox.nominal, fbox.nominal, fbox.scale)
hds = deepcopy(scaled.hd)
nom_s = solve_rhhc(hds, sc, box_nom, s.conic; pce_params = deepcopy(s.pce),
    hpf_optimizer = s.ipopt, nbackground = 0, use_seeds = false)
box0 = ParameterBox(s.box.names, s.box.nominal, s.box.nominal, s.box.nominal, s.box.scale)
nom_l = solve_rhhc(deepcopy(s.hd), s.model, box0, s.conic; pce_params = deepcopy(s.pce),
    hpf_optimizer = s.ipopt, nbackground = 0, use_seeds = false)
gs = Dict(h => I^2 for (h, I) in nom_s["master_solution"])
say("3. √h scaling of the 50 Hz parameters")
say(@sprintf("   nominal: scaled %.6f, physical %.6f (%+.1f %%); scaled solution on the physical model: violation %+.2e at the nominal point",
    nom_s["upper_bound"], nom_l["upper_bound"],
    100(nom_s["upper_bound"] / nom_l["upper_bound"] - 1),
    maximum(violations(hd1, gs, zeros(d))[pos])))
# Family intervals include the eight designs at the fundamental and case harmonics.
# reduce.jl records the subdivision, interior search and containment evidence.
for (n, lo, hi) in zip(fbox.names, fbox.lower, fbox.upper)
    say(@sprintf("   factor %-9s %.4f to %.4f", n, lo, hi))
end
rob_s = solve_rhhc(deepcopy(scaled.hd), sc, fbox, s.conic; pce_params = deepcopy(s.pce),
    hpf_optimizer = s.ipopt, posterior = 256)
pc = posterior_check(hd1, s.model, s.box, rob_s["solution"]; n = 1024)
say(@sprintf("   robust over the factor box: %.6f (%s); on the physical model: worst %+.2e, violated at %.1f %% of 1024 random points (physical robust %.6f)",
    rob_s["lower_bound"], rob_s["signed_off"] ? "signed off" : "NOT signed off",
    pc.max_violation, 100pc.violated_fraction, res["lower_bound"]))
close(io)
