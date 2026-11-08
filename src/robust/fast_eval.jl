################################################################################
# HarmonicPowerModels.jl                                                       #
# Robust Harmonic Hosting Capacity (RHHC): targeted evaluation.                #
#                                                                               #
# The oracle's ascent needs one constraint at many parameter points. The full  #
# evaluation (`moments_at` + `shhc_constraints`) recomputes every constraint   #
# at every harmonic from the network data. Here:                               #
#   - the admittance of the elements the model does not change is assembled   #
#     once per harmonic order and cached; a point only restamps the uncertain #
#     branches (NetworkCache);                                                 #
#   - one constraint needs only its harmonics and one or two rows of the       #
#     transfer impedance, i.e. one LU and one or two solves per harmonic, and #
#     the PCE square of one or two response vectors (ConstraintEvaluator).    #
# The arithmetic follows `shhc_moments`; the values agree to rounding.        #
################################################################################
# Author: Tom Van Acker                                                         #
################################################################################

import LinearAlgebra: lu, transpose

"""
    NetworkCache(hdata, model)

Harmonic admittance at parameter points of `model`, assembled from a cached admittance of the
fixed elements per (harmonic order, sequence) plus the stamps of the uncertain branches.
Thread-safe.
"""
struct NetworkCache{M<:AbstractNetworkModel}
    hdata::Dict
    model::M
    bus_idx::Dict{Int,Int}
    buses::Vector{Int}                         # non-clean bus ids, row order of the reduced Y
    keep::Vector{Int}                          # their matrix indices
    unc::Vector{Int}                           # uncertain branches
    ubr::Dict{Int,Tuple{Int,Int,ComplexF64}}   # branch => (from index, to index, tap)
    yfix::Dict{Tuple{Float64,Symbol},Matrix{ComplexF64}}
    zview::Dict{Int,Dict{String,Any}}          # h => zero-sequence view of the fixed branches
    lk::ReentrantLock
end

_tap(br) = (tm = get(br, "tap", 1.0); tm == 0.0 && (tm = 1.0); tm * cis(get(br, "shift", 0.0)))

function NetworkCache(hdata::Dict, model::AbstractNetworkModel)
    fund    = hdata["nw"]["1"]
    ids     = sort([b["bus_i"] for b in values(fund["bus"])])
    bus_idx = Dict(id => i for (i, id) in enumerate(ids))
    buses   = sort([b["bus_i"] for b in values(fund["bus"]) if !get(b, "is_clean", false)])
    unc     = uncertain_branches(model)
    ubr     = Dict(b => (bus_idx[fund["branch"]["$b"]["f_bus"]], bus_idx[fund["branch"]["$b"]["t_bus"]],
                         _tap(fund["branch"]["$b"])) for b in unc)
    return NetworkCache(hdata, model, bus_idx, buses, [bus_idx[b] for b in buses], unc, ubr,
                        Dict{Tuple{Float64,Symbol},Matrix{ComplexF64}}(),
                        Dict{Int,Dict{String,Any}}(), ReentrantLock())
end

# admittance of the fixed elements at harmonic order h (cached)
function _yfix(nc::NetworkCache, h::Real, sequence::Symbol)
    key = (Float64(h), sequence)
    Y = lock(() -> get(nc.yfix, key, nothing), nc.lk)
    Y === nothing || return Y
    nw  = isinteger(h) && haskey(nc.hdata["nw"], "$(Int(h))") &&
          harmonic_sequence(Int(h)) === sequence ? nc.hdata["nw"]["$(Int(h))"] :
          harmonic_network(nc.hdata, h; sequence)
    ntw = copy(nw)
    ntw["branch"] = Dict(k => v for (k, v) in nw["branch"] if !(parse(Int, k) in nc.unc))
    Y = _harmonic_admittance(ntw, nc.bus_idx)
    lock(() -> (nc.yfix[key] = Y), nc.lk)
    return Y
end

"Full harmonic admittance at the parameter point `p` and the (real) harmonic order `h`."
function fast_admittance(nc::NetworkCache, p, h::Real; sequence::Symbol = :positive)
    Y = copy(_yfix(nc, h, sequence))
    for b in nc.unc
        f, t, tc = nc.ubr[b]
        d = branch_data(nc.model, b, p, h, sequence)
        _stamp_branch!(Y, f, t, d.br_r, d.br_x, complex(d.g_fr, d.b_fr), complex(d.g_to, d.b_to), tc)
    end
    return Y
end

"""
    reduced_admittance(nc::NetworkCache, p, h; sequence = :positive) -> (Y, buses)

As [`reduced_admittance`](@ref) for `(hdata, model, p, h)`, from the cache.
"""
function reduced_admittance(nc::NetworkCache, p, h::Real; sequence::Symbol = :positive)
    Yr = fast_admittance(nc, p, h; sequence)[nc.keep, nc.keep]
    sequence === :zero && (Yr = first(_earth_floating(Yr)))
    return Yr, nc.buses
end

# zero-sequence view of the fixed (not uncertain) branches at the integer harmonic h (cached)
function _zero_branches(nc::NetworkCache, h::Int)
    v = lock(() -> get(nc.zview, h, nothing), nc.lk)
    v === nothing || return v
    nw  = nc.hdata["nw"]["$h"]
    ntw = copy(nw)
    ntw["branch"] = Dict(k => b for (k, b) in nw["branch"] if !(parse(Int, k) in nc.unc))
    v = _zero_sequence_view(ntw)["branch"]
    lock(() -> (nc.zview[h] = v), nc.lk)
    return v
end

_admittance_of(r, x) = abs(r) + abs(x) > 0.0 ? 1.0 / (r + im * x) : zero(ComplexF64)

# terminals (matrix indices, 0 for earth), tap and series admittance of branch b at (p, h);
# at a zero-sequence harmonic a fixed branch is read from its zero-sequence view
function _branch_terminals(ce, b::Int, p, h::Int, bd = nothing)
    f, t, tc = ce.brinfo[b]
    nc = ce.nc
    if harmonic_sequence(h) !== :zero || b in nc.unc
        ys = bd !== nothing && b in nc.unc ? 1.0 / (bd[b].br_r + im * bd[b].br_x) :
             _series_admittance(ce, b, p, h)
        return f, t, tc, ys
    end
    br = _zero_branches(nc, h)["$b"]
    get(br, "earth_fr", false) && (f = 0)
    get(br, "earth_to", false) && (t = 0)
    return f, t, _tap(br), _admittance_of(br["br_r"], br["br_x"])
end

"""
    ConstraintEvaluator(hdata, model, box)

Normalized violation of one chance constraint at a latent point, computed from only the
harmonics and transfer-impedance rows that constraint needs. `hdata` must carry the PCE data.
"""
struct ConstraintEvaluator{M<:AbstractNetworkModel}
    nc::NetworkCache{M}
    box::ParameterBox
    pb::Any
    K::Vector{Int}
    psi_norms::Vector{Float64}
    E::Dict{Int,Matrix{ComplexF64}}            # h => injected-current PCE rows w_d(cos + j sin)
    upos::Dict{Int,Vector{Int}}                # h => position in the reduced Y of each row's bus
    pos::Dict{Int,Int}                         # matrix index => position in the reduced Y
    cons::Vector{ShhcConstraint}               # nominal records (kind, key, harmonics, limit, λ)
    brinfo::Dict{Int,Tuple{Int,Int,ComplexF64}}   # branch => (from, to, tap) for :rmsi
    deps::Vector{Vector{Int}}                  # coordinate => uncertain branches it changes
    recomb::Bool                               # false: worst-case alignment (no recombination)
    W::Dict{Int,Vector{Float64}}               # h => allocation weight w_d of each unit row
end

function ConstraintEvaluator(hdata::Dict, model::AbstractNetworkModel, box::ParameterBox)
    nc   = NetworkCache(hdata, model)
    pce  = hdata["pce"]
    m0   = moments_at(hdata, model, nominal_point(box))
    pos  = Dict(i => j for (j, i) in enumerate(nc.keep))
    wts  = fairness_weights(hdata, pce["harmonic_units"])
    E    = Dict{Int,Matrix{ComplexF64}}(); upos = Dict{Int,Vector{Int}}()
    recomb = get(pce, "recombination", true)
    hunits = Set(pce["harmonic_units"]); W = Dict{Int,Vector{Float64}}()
    for h in m0.harm_nw_ints
        haskey(pce["cos_phi"], h) || continue
        rows = Vector{ComplexF64}[]; ps = Int[]; ws = Float64[]
        for (d_str, load) in hdata["nw"]["$h"]["load"]
            d = parse(Int, d_str); bd = nc.bus_idx[load["load_bus"]]
            haskey(pos, bd) || continue                      # unit at an ideal ground
            w = get(wts, d, 0.0)
            # without recombination only the harmonic units count (normal loads inject nothing)
            recomb || (d in hunits && w != 0.0) || continue
            push!(rows, w .* (pce["cos_phi"][h][d] .+ im .* pce["sin_phi"][h][d])); push!(ps, pos[bd])
            push!(ws, w)
        end
        E[h] = isempty(rows) ? zeros(ComplexF64, 0, length(m0.K)) : permutedims(reduce(hcat, rows))
        upos[h] = ps; W[h] = ws
    end
    fund   = hdata["nw"]["1"]["branch"]
    brinfo = Dict(parse(Int, k) => (nc.bus_idx[b["f_bus"]], nc.bus_idx[b["t_bus"]], _tap(b))
                  for (k, b) in fund)
    deps   = _dependencies(nc, box, minimum(m0.harm_nw_ints))
    return ConstraintEvaluator(nc, box, m0.pb, collect(m0.K), collect(Float64, m0.psi_norms), E,
                               upos, pos, shhc_constraints(hdata, m0), brinfo, deps, recomb, W)
end

# PCE coefficients (re, im) of the voltage response at matrix index n, from the factorised
# transpose of the reduced Y (row n of its inverse = transfer impedances from the units)
function _response(ce::ConstraintEvaluator, Ft, n::Int, h::Int)
    nk = length(ce.K)
    haskey(ce.pos, n) || return zeros(nk), zeros(nk)    # ideal ground: V = 0
    e = zeros(ComplexF64, length(ce.pos)); e[ce.pos[n]] = 1.0
    x = Ft \ e
    C = isempty(ce.upos[h]) ? zeros(ComplexF64, nk) : vec(transpose(x[ce.upos[h]]) * ce.E[h])
    return real.(C), imag.(C)
end

# row n of the inverse of the reduced Y (zero for an ideal ground or for earth, n = 0)
function _zrow(ce::ConstraintEvaluator, Ft, n::Int)
    haskey(ce.pos, n) || return zeros(ComplexF64, length(ce.pos))
    e = zeros(ComplexF64, length(ce.pos)); e[ce.pos[n]] = 1.0
    return Ft \ e
end

# worst-case (no recombination) coefficient: (Σ_d w_d |u_d|)² as the mean, zero variance
function _worst_case_coef(ce::ConstraintEvaluator, u::AbstractVector, h::Int)
    c = zeros(length(ce.K)); c[1] = sum(ce.W[h] .* abs.(u); init = 0.0)^2
    return c
end

# series admittance of branch b at (p, h): from the model if uncertain, else from the data
function _series_admittance(ce::ConstraintEvaluator, b::Int, p, h::Int)
    nc = ce.nc
    if b in nc.unc
        d = branch_data(nc.model, b, p, h, harmonic_sequence(h)); r, x = d.br_r, d.br_x
    else
        br = nc.hdata["nw"]["$h"]["branch"]["$b"]; r, x = br["br_r"], br["br_x"]
    end
    return 1.0 / (r + im * x)
end

"""
    fast_violation(ce, j, ξ, g) -> Float64

Normalized violation of constraint `j` (index in `shhc_constraints`) at the latent point `ξ`
for the hosting variables `g` (h ⇒ g_h = I_h²); equals [`normalized_violation`](@ref) of the
full evaluation to rounding.
"""
function fast_violation(ce::ConstraintEvaluator, j::Int, ξ, g)
    c  = ce.cons[j]
    p  = latent(ce.box, ξ)
    coef = Dict{Int,Vector{Float64}}()
    for h in keys(c.coef)
        Yr = reduced_admittance(ce.nc, p, h; sequence = harmonic_sequence(h))[1]
        Ft = lu(transpose(Yr))
        if !ce.recomb
            if c.kind === :rmsi
                f, t, tc, ys = _branch_terminals(ce, c.key, p, h)
                u = ys .* ((1.0 / tc) .* _zrow(ce, Ft, f) .- _zrow(ce, Ft, t))
            else
                bus = c.kind === :ihd ? c.key[2] : c.key
                u = _zrow(ce, Ft, ce.nc.bus_idx[bus])
            end
            coef[h] = _worst_case_coef(ce, u[ce.upos[h]], h)
            continue
        end
        if c.kind === :rmsi
            f, t, tc, ys = _branch_terminals(ce, c.key, p, h); gy, by = real(ys), imag(ys)
            arf, aif = _response(ce, Ft, f, h); art, ait = _response(ce, Ft, t, h)
            ct = 1.0 / tc; cr, ci = real(ct), imag(ct)
            dVr = (cr .* arf .- ci .* aif) .- art
            dVi = (cr .* aif .+ ci .* arf) .- ait
            Ibr = gy .* dVr .- by .* dVi
            Ibi = gy .* dVi .+ by .* dVr
            coef[h] = galerkin_product(ce.pb, Ibr, Ibr) .+ galerkin_product(ce.pb, Ibi, Ibi)
        else
            bus = c.kind === :ihd ? c.key[2] : c.key
            ar, ai = _response(ce, Ft, ce.nc.bus_idx[bus], h)
            coef[h] = galerkin_product(ce.pb, ar, ar) .+ galerkin_product(ce.pb, ai, ai)
        end
    end
    F = ShhcConstraint(c.kind, c.key, coef, c.limit, c.λ)
    return normalized_violation(F, g, (K = ce.K, psi_norms = ce.psi_norms))
end

# ──────────────────────────────────────────────────────────────────────────────
# Analytic gradient (adjoint)
# ──────────────────────────────────────────────────────────────────────────────
#
# For a constraint whose response at harmonic h is C = Σ_d z_d E_d, with z the row of
# Z = Y⁻¹ of its bus (or, for a branch current, y_s·(t⁻¹·z_f − z_t)), the derivative of
# F along a parameter follows from dZ = −Z dY Z:
#
#     dF = Σ_h −Re(rᵀ (∂Y/∂p) v) [+ Re(dy_s · Σ_k conj(G_k) D_k) for an uncertain branch]
#
# with r the response row, G = ∂F/∂C (complex, from the PCE square), q = E·conj(G) and the
# adjoint v = Y⁻¹ s, s the scatter of q onto the unit buses. One LU and two solves per
# harmonic give the whole gradient; ∂Y/∂p only involves the stamps of the uncertain
# branches the parameter changes, whose data are differentiated by central differences
# or ForwardDiff dual numbers (no linear algebra).

# ∂(z_k)/∂x of z = galerkin_product(pb, x, x), contracted with β: Σ_k β_k ∂z_k/∂x
function _galerkin_grad(pb, β::AbstractVector, x::AbstractVector)
    gr = zeros(length(x))
    for (j, l, k, v) in pb.M_nz
        if j == l
            gr[j] += 2v * β[k] * x[j]
        else
            gr[j] += 2v * β[k] * x[l]
            gr[l] += 2v * β[k] * x[j]
        end
    end
    return gr
end

# uncertain branches whose data change with each latent coordinate (probed at the nominal point)
function _dependencies(nc::NetworkCache, box::ParameterBox, h::Int)
    d   = length(box)
    deps = [Int[] for _ in 1:d]
    p0  = latent(box, zeros(d))
    for i in 1:d
        e = zeros(d); e[i] = 0.5
        p₊, p₋ = latent(box, e), latent(box, -e)
        for b in nc.unc
            seq = harmonic_sequence(h)
            branch_data(nc.model, b, p₊, h, seq) == branch_data(nc.model, b, p₋, h, seq) &&
                branch_data(nc.model, b, p₊, h, seq) == branch_data(nc.model, b, p0, h, seq) ||
                push!(deps[i], b)
        end
    end
    return deps
end

# HPM tag and maximum number of latent partials in one branch evaluation
struct _BranchDerivativeTag end
const _BRANCH_DERIVATIVE_CHUNK = 32

function _check_branch_derivative(branch_derivative)
    (branch_derivative === :central || branch_derivative === :dual) || throw(ArgumentError(
        "branch_derivative must be :central or :dual, got $(repr(branch_derivative))"))
    return nothing
end

# preserve the original diagnostic when a model cannot evaluate a dual point
function _dual_branch_call(f, b)
    try
        return f()
    catch err
        err isa MethodError || rethrow()
        error("branch $b: the data do not accept dual numbers; use " *
              "branch_derivative = :central. Original error: " * sprint(showerror, err))
    end
end

function _dual_branch_data(model::AbstractNetworkModel, branches, p, h, sequence)
    return Dict(b => _dual_branch_call(() -> branch_data(model, b, p, h, sequence), b)
                for b in branches)
end

function _dual_branch_data(model::ExactPiModel, branches, p, h, sequence)
    groups = Vector{Int}[]
    for b in branches
        k = findfirst(bs -> model.branches[first(bs)].zy === model.branches[b].zy, groups)
        k === nothing ? push!(groups, [b]) : push!(groups[k], b)
    end
    data = Dict{Int,NamedTuple}()
    for bs in groups
        b = first(bs)
        z, y = _dual_branch_call(() -> model.branches[b].zy(p, model.f1 * h, sequence), b)
        for b in bs
            br = model.branches[b]
            data[b] = _dual_branch_call(b) do
                len = br.len isa Function ? br.len(p) : br.len
                _exact_pi(br, len, z, y)
            end
        end
    end
    return data
end

# only the partials enter the adjoint; its values remain those of the plain forward pass
_branch_partial(x::Real, i) = 0.0
_branch_partial(x::ForwardDiff.Dual, i) = ForwardDiff.partials(x)[i]

function _dual_branch_derivatives(ce, ξ, h, sequence)
    d = length(ξ)
    data = Dict{Tuple{Int,Int},NamedTuple}()
    for lo in 1:_BRANCH_DERIVATIVE_CHUNK:d
        coords = lo:min(lo + _BRANCH_DERIVATIVE_CHUNK - 1, d)
        branches = sort(unique(vcat(ce.deps[coords]...)))
        isempty(branches) && continue
        n = length(coords)
        ξd = [ForwardDiff.Dual{_BranchDerivativeTag}(ξ[i],
                  ntuple(k -> i == coords[k] ? 1.0 : 0.0, n)) for i in 1:d]
        bd = _dual_branch_data(ce.nc.model, branches, latent(ce.box, ξd), h, sequence)
        for (k, i) in enumerate(coords), b in ce.deps[i]
            data[(i, b)] = map(x -> _branch_partial(x, k), bd[b])
        end
    end
    return data
end

"""
    fast_violation_gradient(ce, j, ξ, g; δ = 1e-6, branch_derivative = :central) -> (value, gradient)

Normalized violation of constraint `j` at the latent point `ξ` and its gradient with respect
to `ξ`, by the adjoint of the transfer impedances (see the comment above). The branch data of
the uncertain branches are differentiated by central differences with step `δ`
(`branch_derivative = :central`) or ForwardDiff dual numbers (`:dual`).

`:dual` needs branch data that carry the partials of a dual point through. A model that drops
them without an error, for example by solving an equation by bisection on plain bounds,
indexing a table by a value, or taking the value of a dual number, gives a wrong gradient.
Compare the two options on the model before relying on `:dual`.
"""
function fast_violation_gradient(ce::ConstraintEvaluator, j::Int, ξ::AbstractVector, g;
                                 δ::Real = 1e-6, branch_derivative = :central)
    _check_branch_derivative(branch_derivative)
    c   = ce.cons[j]
    nc  = ce.nc
    d   = length(ξ)
    p   = latent(ce.box, ξ)
    hs  = sort(collect(keys(c.coef)))
    nk  = length(ce.K)
    deps = ce.deps

    # forward pass per harmonic: response row r, PCE response C and, for a branch, D
    fw = Dict{Int,Any}()
    coef = Dict{Int,Vector{Float64}}()
    for h in hs
        seq = harmonic_sequence(h)
        bd  = Dict(b => branch_data(nc.model, b, p, h, seq) for b in nc.unc)
        Y   = copy(_yfix(nc, h, seq))
        for b in nc.unc
            f, t, tc = nc.ubr[b]; x = bd[b]
            _stamp_branch!(Y, f, t, x.br_r, x.br_x, complex(x.g_fr, x.b_fr), complex(x.g_to, x.b_to), tc)
        end
        Yr  = Y[nc.keep, nc.keep]
        seq === :zero && (Yr = first(_earth_floating(Yr)))
        F   = lu(Yr)
        row(n) = haskey(ce.pos, n) ? (e = zeros(ComplexF64, length(ce.pos)); e[ce.pos[n]] = 1.0;
                                       transpose(F) \ e) : zeros(ComplexF64, length(ce.pos))
        resp(rv) = isempty(ce.upos[h]) ? zeros(ComplexF64, nk) : vec(transpose(rv[ce.upos[h]]) * ce.E[h])
        if !ce.recomb
            # worst-case alignment: u_d per unit, the sum of magnitudes S_h, D_d = t⁻¹z_f − z_t
            if c.kind === :rmsi
                f, t, tc, ys = _branch_terminals(ce, c.key, p, h, bd)
                Dv = (1.0 / tc) .* row(f) .- row(t)
                r  = ys .* Dv
                Du = Dv[ce.upos[h]]
            else
                bus = c.kind === :ihd ? c.key[2] : c.key
                r   = row(nc.bus_idx[bus]); Du = nothing
            end
            u  = r[ce.upos[h]]
            Sh = sum(ce.W[h] .* abs.(u); init = 0.0)
            fw[h] = (F = F, r = r, C = nothing, D = Du, ys = nothing, bd = bd, u = u, Sh = Sh)
            coef[h] = _worst_case_coef(ce, u, h)
            continue
        end
        if c.kind === :rmsi
            f, t, tc, ys = _branch_terminals(ce, c.key, p, h, bd)
            ct = 1.0 / tc
            xf, xt = row(f), row(t)
            r  = ys .* (ct .* xf .- xt)
            D  = resp(ct .* xf .- xt)
            fw[h] = (F = F, r = r, C = ys .* D, D = D, ys = ys, bd = bd)
        else
            bus = c.kind === :ihd ? c.key[2] : c.key
            r   = row(nc.bus_idx[bus])
            fw[h] = (F = F, r = r, C = resp(r), D = nothing, ys = nothing, bd = bd)
        end
        C = fw[h].C
        coef[h] = galerkin_product(ce.pb, real.(C), real.(C)) .+ galerkin_product(ce.pb, imag.(C), imag.(C))
    end

    # value and ∂F/∂α_h
    ES  = sum(coef[h][1] * g[h] for h in hs)
    Kh  = [k for k in ce.K if k != 0]
    S   = Dict(k => sum(coef[h][k+1] * g[h] for h in hs) for k in Kh)
    var = sum(ce.psi_norms[k+1] * S[k]^2 for k in Kh; init = 0.0)
    sd  = sqrt(var)
    Fv  = c.λ <= 0.0 ? ES : ES + c.λ * sd
    L   = c.limit
    val = L > 0.0 ? Fv / L - 1.0 : Fv - L
    scale = L > 0.0 ? 1.0 / L : 1.0

    grad = zeros(d)
    for h in hs
        β = zeros(nk); β[1] = g[h]
        if c.λ > 0.0 && sd > 0.0
            for k in Kh
                β[k+1] = c.λ * ce.psi_norms[k+1] * S[k] * g[h] / sd
            end
        end
        if ce.recomb
            C  = fw[h].C
            G  = _galerkin_grad(ce.pb, β, real.(C)) .+ im .* _galerkin_grad(ce.pb, β, imag.(C))
            q  = ce.E[h] * conj.(G)
        else
            # ∂(g_h S_h²)/∂u_d contracted as Re(c_d du_d), c_d = 2 g_h S_h w_d conj(u_d)/|u_d|
            u  = fw[h].u
            q  = [abs(u[a]) > 0.0 ? 2 * g[h] * fw[h].Sh * ce.W[h][a] * conj(u[a]) / abs(u[a]) :
                  zero(ComplexF64) for a in eachindex(u)]
        end
        s  = zeros(ComplexF64, length(ce.pos))
        for (a, pa) in enumerate(ce.upos[h]); s[pa] += q[a]; end
        v  = fw[h].F \ s
        r  = fw[h].r
        seq = harmonic_sequence(h)
        db = branch_derivative === :dual ? _dual_branch_derivatives(ce, ξ, h, seq) : nothing
        for i in 1:d, b in deps[i]
            x  = fw[h].bd[b]
            ys = 1.0 / (x.br_r + im * x.br_x)
            if branch_derivative === :central
                e = zeros(d); e[i] = δ
                x₊ = branch_data(nc.model, b, latent(ce.box, ξ .+ e), h, seq)
                x₋ = branch_data(nc.model, b, latent(ce.box, ξ .- e), h, seq)
                dys  = -ys^2 * ((x₊.br_r - x₋.br_r) + im * (x₊.br_x - x₋.br_x)) / (2δ)
                dyfr = ((x₊.g_fr - x₋.g_fr) + im * (x₊.b_fr - x₋.b_fr)) / (2δ)
                dyto = ((x₊.g_to - x₋.g_to) + im * (x₊.b_to - x₋.b_to)) / (2δ)
            else
                dx = db[(i, b)]
                dys  = -ys^2 * (dx.br_r + im * dx.br_x)
                dyfr = dx.g_fr + im * dx.b_fr
                dyto = dx.g_to + im * dx.b_to
            end
            f, t, tc = nc.ubr[b]
            # rᵀ dY v over the four stamp entries, kept buses only
            acc = zero(ComplexF64)
            for (a, bb, dy) in ((f, f, (dys + dyfr) / abs2(tc)), (t, t, dys + dyto),
                                (f, t, -dys / conj(tc)), (t, f, -dys / tc))
                (haskey(ce.pos, a) && haskey(ce.pos, bb)) || continue
                acc += r[ce.pos[a]] * dy * v[ce.pos[bb]]
            end
            grad[i] -= real(acc)
            # a branch-current constraint on an uncertain branch: its own y_s moves too
            if c.kind === :rmsi && b == c.key
                grad[i] += ce.recomb ? real(dys * sum(conj.(G) .* fw[h].D)) :
                                       real(dys * sum(q .* fw[h].D))
            end
        end
    end
    return val, scale .* grad
end
