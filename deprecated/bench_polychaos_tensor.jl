# bench_polychaos_tensor.jl
# Pitfall quantification for the PolyChaos.jl migration (scratch script).
#   A. numerical agreement native vs PolyChaos at industrial size (N=4)
#   B. default-quadrature exactness of Tensor(3, op)
#   C. scale test at IEEE-30 size (N=21, P=253): build time, memory, access
# Run:  julia --project bench_polychaos_tensor.jl [--skip-large]

using PolyChaos, LinearAlgebra, Printf

include(joinpath(@__DIR__, "src/core/pce.jl"))   # native implementation

const α̂ = 4.0414
const β̂ = 4.7762
const δ = 2
const a = β̂ - 1.0
const b = α̂ - 1.0

function native_tensor(N)
    mi      = generate_multi_indices(N, δ)
    nodes, wts = gauss_jacobi_nodes_weights(3δ + 2, a, b)
    norms_u = [beta_poly_norm_sq(n, a, b) for n in 0:δ]
    t = @elapsed M = compute_M_tensor(mi, norms_u, nodes, wts, a, b)
    return mi, M, t
end

pc_mop(N; Nrec = 3δ + 2) =
    MultiOrthoPoly([JacobiOrthoPoly(δ, a, b; Nrec = Nrec) for _ in 1:N], δ)

# computeSP2(::MultiOrthoPoly) is broken in PolyChaos 1.1 (accesses op.β);
# build the multivariate norms from the univariate ones and mop.ind instead.
function sp2_multi(mop)
    sp2u = computeSP2(mop.uni[1])
    return [prod(sp2u[mop.ind[r, u] + 1] for u in 1:size(mop.ind, 2))
            for r in 1:size(mop.ind, 1)]
end

# map native multi-index ordering -> PolyChaos row index
function index_map(mi_native, ind_pc)
    d = Dict{Vector{Int},Int}()
    for r in 1:size(ind_pc, 1)
        d[vec(ind_pc[r, :])] = r
    end
    return [d[m] for m in mi_native]
end

# ── A. correctness, industrial size ─────────────────────────────────────────
println("=== A. correctness  (N=4, industrial size) ===")
mi4, M4, _ = native_tensor(4)
P4   = length(mi4)
mop4 = pc_mop(4)
println("P native = $P4, P polychaos = $(mop4.dim)")
sp2  = sp2_multi(mop4)
t3   = Tensor(3, mop4)
perm = index_map(mi4, mop4.ind)
maxdiff = 0.0
for j in 1:P4, l in 1:P4, k in 1:P4
    pc = t3.get([perm[j]-1, perm[l]-1, perm[k]-1]) /
         sqrt(sp2[perm[j]] * sp2[perm[l]] * sp2[perm[k]])
    global maxdiff = max(maxdiff, abs(pc - M4[j, l, k]))
end
@printf "max |M_pc - M_native| over all %d^3 entries = %.3e\n\n" P4 maxdiff

# ── B. default-quadrature exactness ─────────────────────────────────────────
println("=== B. Tensor exactness vs Nrec  (univariate, δ=$δ) ===")
for nrec in [δ + 1, 2δ + 1, 3δ + 2, 20]
    op = JacobiOrthoPoly(δ, a, b; Nrec = nrec)
    try
        t3u = Tensor(3, op)
        @printf "Nrec = %2d  →  <ψ2ψ2ψ2> = %.12f   (quad nodes: %d)\n" nrec t3u.get([2, 2, 2]) length(op.quad.nodes)
    catch e
        @printf "Nrec = %2d  →  THROWS: %s   (quad nodes: %d)\n" nrec sprint(showerror, e) length(op.quad.nodes)
    end
end
println()

# ── C. scale test, IEEE-30 size ─────────────────────────────────────────────
function run_large()
    println("=== C. scale  (N=21, IEEE-30 size) ===")
    N = 21
    mi21 = generate_multi_indices(N, δ)
    P = length(mi21)
    @printf "P = %d,  dense tensor memory = %.1f MB\n" P P^3 * 8 / 1024^2

    nodes, wts = gauss_jacobi_nodes_weights(3δ + 2, a, b)
    norms_u    = [beta_poly_norm_sq(n, a, b) for n in 0:δ]
    t_nat = @elapsed M21 = compute_M_tensor(mi21, norms_u, nodes, wts, a, b)
    nnz_M = count(x -> abs(x) > 1e-14, M21)
    @printf "native compute_M_tensor      : %7.2f s\n" t_nat
    @printf "nonzeros                     : %d / %d  (%.3f %%)\n" nnz_M P^3 100nnz_M / P^3

    t_set = @elapsed mop21 = pc_mop(N)
    t_sp2 = @elapsed sp221 = sp2_multi(mop21)
    t_t3  = @elapsed t321  = Tensor(3, mop21)
    @printf "PolyChaos MultiOrthoPoly     : %7.2f s\n" t_set
    @printf "PolyChaos computeSP2         : %7.2f s\n" t_sp2
    @printf "PolyChaos Tensor(3, mop)     : %7.2f s\n" t_t3
    @printf "Tensor storage (summarysize) : %7.1f MB\n" Base.summarysize(t321) / 1024^2

    # spot-check agreement on 200 random triples
    perm21 = index_map(mi21, mop21.ind)
    md = 0.0
    for _ in 1:200
        j, l, k = rand(1:P), rand(1:P), rand(1:P)
        pc = t321.get([perm21[j]-1, perm21[l]-1, perm21[k]-1]) /
             sqrt(sp221[perm21[j]] * sp221[perm21[l]] * sp221[perm21[k]])
        md = max(md, abs(pc - M21[j, l, k]))
    end
    @printf "spot-check max diff (200 random triples): %.3e\n" md

    # access cost: full P^3 sweep through t.get (what a naive port would do)
    t_get = @elapsed begin
        s = 0.0
        idx = Vector{Int}(undef, 3)
        for j in 0:P-1, l in 0:P-1, k in 0:P-1
            idx[1] = j; idx[2] = l; idx[3] = k
            s += t321.get(idx)
        end
        s
    end
    @printf "full P^3 sweep via t.get     : %7.2f s\n" t_get

    # Galerkin square (compute_alpha inner kernel), dense, one bus
    vr, vi = randn(P), randn(P)
    alpha  = zeros(P)
    t_alpha = @elapsed begin
        for k in 1:P, k1 in 1:P, k2 in 1:P
            alpha[k] += M21[k1, k2, k] * (vr[k1] * vr[k2] + vi[k1] * vi[k2])
        end
    end
    @printf "dense Galerkin square (1 bus): %7.2f s\n" t_alpha

    # same kernel restricted to the nonzero pattern (sparse-style iteration)
    nz = [(k1, k2, k, M21[k1, k2, k]) for k in 1:P for k1 in 1:P for k2 in 1:P
          if abs(M21[k1, k2, k]) > 1e-14]
    alpha2 = zeros(P)
    t_alpha_sp = @elapsed begin
        for (k1, k2, k, v) in nz
            alpha2[k] += v * (vr[k1] * vr[k2] + vi[k1] * vi[k2])
        end
    end
    @printf "sparse Galerkin square (1 bus): %6.4f s   (max dev %.2e)\n" t_alpha_sp maximum(abs.(alpha2 - alpha))
end

"--skip-large" in ARGS || run_large()
