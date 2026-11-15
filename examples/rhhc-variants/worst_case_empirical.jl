################################################################################
# examples/rhhc-variants/worst_case_empirical.jl                               #
# A-posteriori check of the Gaussian factor at the worst cases of the robust   #
# solution on `ieee30-full-exact` (23 coordinates, harmonics 5–25, exact       #
# moments, ε = 0.05). For every constraint with a scenario set beyond the      #
# nominal point, and every point of that set, at the robust hosting capacity:  #
#   - the empirical factor λ* = (q₁₋ε − E)/σ of the sampled output,            #
#   - the sampled (1−ε) quantile relative to the limit, q/a, and               #
#   - the violation probability P(X > a), by bisection on ε with the same      #
#     samples (fixed seed), so that q at that probability equals the limit.    #
#                                                                               #
# SHHC_HMAX=25 julia -t 10 --project=examples/rhhc-variants \                   #
#     examples/rhhc-variants/worst_case_empirical.jl                            #
# Output: output/ieee30-full-exact-worst-case-empirical.txt                     #
################################################################################

isdefined(@__MODULE__, :STUDIES) || include(joinpath(@__DIR__, "builders.jl"))
using Printf, Statistics, Distributions

const OUT = joinpath(STUDIES, "output", "ieee30-full-exact-worst-case-empirical.txt")
io = open(OUT, "w")
say(args...) = (println(io, args...); println(args...); flush(io))

s = study("cases/ieee30-full-exact.toml")
header(io, s)
hd  = deepcopy(s.hd)                       # prepared in place by solve_rhhc
res = solve_rhhc(hd, s.model, s.box, s.conic; pce_params = deepcopy(s.pce),
                 hpf_optimizer = s.ipopt, posterior = 1024)
g    = Dict(h => I^2 for (h, I) in res["solution"])
ε    = hd["pce"]["epsilon"]
λg   = quantile(Normal(), 1 - ε)
dict = res["dictionary"]
say(@sprintf("robust %.6f, %d iterations, %s; Gaussian factor %.4f at ε = %.2f",
    res["lower_bound"], res["iterations"], res["signed_off"] ? "signed off" : "NOT signed off", λg, ε))

group(kind) = kind in (:thd, :rmsv) ? :agg : kind
const NSAMP = 20_000

# sampled quantile of the constraint of `kind`/`key` at the point hdata `hp`, at risk level e
function sampled(hp, kind, key, e)
    hp["pce"]["epsilon"] = e
    λs = HPM.empirical_lambdas(hp, g; nsamp = NSAMP, seed = 1)
    # the IHD quantile is per unit of g_h (|V|²/I²); the aggregates include g already
    scale = kind === :ihd ? g[key[1]] : 1.0
    return scale * λs[:q][group(kind)][key], λs[group(kind)][key]
end

rows = Tuple{Any,Int,Float64,Float64,Float64,Float64}[]   # key, point, F/a, λ*, q/a, P(X>a)
for (key, idx) in sort(collect(res["scenarios"]); by = string)
    kind, ckey = key
    for k in idx
        k == 1 && continue                                     # the nominal point
        hp = HPM.hdata_at(hd, s.model, latent(s.box, dict.points[k]))
        m  = HPM.shhc_moments(hp)
        c  = only(filter(c -> c.kind == kind && c.key == ckey, HPM.shhc_constraints(hp, m)))
        Fa = HPM.constraint_value(c, g, m) / c.limit
        q, λstar = sampled(hp, kind, ckey, ε)
        # violation probability: q₁₋ₑ(X) = a, bisection on e with the same samples
        lo, hi = 1e-5, 0.5
        p = if sampled(hp, kind, ckey, lo)[1] <= c.limit
            0.0
        else
            for _ in 1:14
                mid = sqrt(lo * hi)
                sampled(hp, kind, ckey, mid)[1] > c.limit ? (lo = mid) : (hi = mid)
            end
            sqrt(lo * hi)
        end
        push!(rows, (key, k, Fa, λstar, q / c.limit, p))
        say(@sprintf("   %-22s point %5d (%-10s)  F/a %.4f  λ* %.4f  q/a %.4f  P(X > a) %.4f",
            string(key), k, dict.labels[k], Fa, λstar, q / c.limit, p))
    end
end
binding = filter(r -> r[3] > 0.99, rows)
say(@sprintf("%d worst-case points over %d constraints; %d with the Gaussian bound within 1 %% of the limit",
    length(rows), length(res["scenarios"]), length(binding)))
if !isempty(rows)
    say(@sprintf("all points:     λ* %.3f to %.3f, q/a %.4f to %.4f, P(X > a) %.4f to %.4f (ε = %.2f)",
        extrema(r[4] for r in rows)..., extrema(r[5] for r in rows)..., extrema(r[6] for r in rows)..., ε))
end
if !isempty(binding)
    say(@sprintf("binding points: λ* %.3f to %.3f, q/a %.4f to %.4f, P(X > a) %.4f to %.4f",
        extrema(r[4] for r in binding)..., extrema(r[5] for r in binding)..., extrema(r[6] for r in binding)...))
end
close(io)
