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
                        Dict{Tuple{Float64,Symbol},Matrix{ComplexF64}}(), ReentrantLock())
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
    return fast_admittance(nc, p, h; sequence)[nc.keep, nc.keep], nc.buses
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
end

function ConstraintEvaluator(hdata::Dict, model::AbstractNetworkModel, box::ParameterBox)
    nc   = NetworkCache(hdata, model)
    pce  = hdata["pce"]
    m0   = moments_at(hdata, model, nominal_point(box))
    pos  = Dict(i => j for (j, i) in enumerate(nc.keep))
    wts  = fairness_weights(hdata, pce["harmonic_units"])
    E    = Dict{Int,Matrix{ComplexF64}}(); upos = Dict{Int,Vector{Int}}()
    for h in m0.harm_nw_ints
        haskey(pce["cos_phi"], h) || continue
        rows = Vector{ComplexF64}[]; ps = Int[]
        for (d_str, load) in hdata["nw"]["$h"]["load"]
            d = parse(Int, d_str); bd = nc.bus_idx[load["load_bus"]]
            haskey(pos, bd) || continue                      # unit at an ideal ground
            w = get(wts, d, 0.0)
            push!(rows, w .* (pce["cos_phi"][h][d] .+ im .* pce["sin_phi"][h][d])); push!(ps, pos[bd])
        end
        E[h] = isempty(rows) ? zeros(ComplexF64, 0, length(m0.K)) : permutedims(reduce(hcat, rows))
        upos[h] = ps
    end
    fund   = hdata["nw"]["1"]["branch"]
    brinfo = Dict(parse(Int, k) => (nc.bus_idx[b["f_bus"]], nc.bus_idx[b["t_bus"]], _tap(b))
                  for (k, b) in fund)
    return ConstraintEvaluator(nc, box, m0.pb, collect(m0.K), collect(Float64, m0.psi_norms), E,
                               upos, pos, shhc_constraints(hdata, m0), brinfo)
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
        if c.kind === :rmsi
            f, t, tc = ce.brinfo[c.key]
            ys = _series_admittance(ce, c.key, p, h); gy, by = real(ys), imag(ys)
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
