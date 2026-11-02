################################################################################
# HarmonicPowerModels.jl                                                       #
# Robust Harmonic Hosting Capacity (RHHC): parametric network model.           #
#                                                                               #
# A network model maps a point p of a box of uncertain parameters to the       #
# per-unit π-model data of the uncertain branches at any harmonic order h      #
# (real, so also between harmonics). From it follow the harmonic networks at p #
# (input of `shhc_moments`) and the harmonic admittance at p and h (input of   #
# the modal analysis that samples the box).                                    #
################################################################################
# Author: Tom Van Acker                                                         #
################################################################################

# ──────────────────────────────────────────────────────────────────────────────
# Parameter box
# ──────────────────────────────────────────────────────────────────────────────

"""
    ParameterBox(entries...)

Box of uncertain network parameters. Each entry is `name => (nominal, lower, upper)` or
`name => (nominal, lower, upper, scale)`, with `lower ≤ nominal ≤ upper` and `scale` either
`:linear` (default) or `:log` (affine in `log(value)`, for parameters spanning orders of
magnitude, e.g., earth resistivity). The latent coordinates ξ ∈ [-1, 1]^d map onto the box
through [`latent`](@ref).
"""
struct ParameterBox
    names::Vector{Symbol}
    nominal::Vector{Float64}
    lower::Vector{Float64}
    upper::Vector{Float64}
    scale::Vector{Symbol}
    function ParameterBox(names, nominal, lower, upper, scale)
        allunique(names) || throw(ArgumentError("parameter names must be unique: $names"))
        for (n, x, lo, hi, sc) in zip(names, nominal, lower, upper, scale)
            all(isfinite, (x, lo, hi)) ||
                throw(DomainError((lo, x, hi), "$n: bounds must be finite"))
            lo <= x <= hi ||
                throw(DomainError((lo, x, hi), "$n: need lower ≤ nominal ≤ upper"))
            sc in (:linear, :log) ||
                throw(ArgumentError("$n: scale must be :linear or :log, got :$sc"))
            sc === :log && lo <= 0.0 &&
                throw(DomainError(lo, "$n: a :log scale needs a positive lower bound"))
        end
        return new(collect(Symbol, names), collect(Float64, nominal), collect(Float64, lower),
                   collect(Float64, upper), collect(Symbol, scale))
    end
end

function ParameterBox(entries::Pair{Symbol}...)
    names = first.(entries)
    vals  = last.(entries)
    return ParameterBox(names, getindex.(vals, 1), getindex.(vals, 2), getindex.(vals, 3),
                        [length(v) ≥ 4 ? v[4] : :linear for v in vals])
end

Base.length(box::ParameterBox) = length(box.names)

_named(box::ParameterBox, values) = NamedTuple{Tuple(box.names)}(Tuple(values))

"Nominal point of `box` as a named tuple."
nominal_point(box::ParameterBox) = _named(box, box.nominal)

"""
    latent(box, ξ) -> NamedTuple

Map the latent point `ξ ∈ [-1, 1]^d` onto `box`: affine on a `:linear` scale, affine in
`log` on a `:log` scale. Generic in the element type of `ξ` (e.g., dual numbers).
"""
function latent(box::ParameterBox, ξ)
    length(ξ) == length(box) ||
        throw(DimensionMismatch("ξ has length $(length(ξ)), the box has $(length(box))"))
    # convex combination ((1-ξ)·lo + (1+ξ)·hi)/2: exact at both bounds, so vertices
    # evaluate the box's own bounds rather than a rounded neighbor
    vals = map(1:length(box)) do i
        lo, hi = box.lower[i], box.upper[i]
        if box.scale[i] === :linear
            ((1 - ξ[i]) * lo + (1 + ξ[i]) * hi) / 2
        else
            exp(((1 - ξ[i]) * log(lo) + (1 + ξ[i]) * log(hi)) / 2)
        end
    end
    return _named(box, vals)
end

# ──────────────────────────────────────────────────────────────────────────────
# Network models
# ──────────────────────────────────────────────────────────────────────────────

"""
    AbstractNetworkModel

A parametric model of the uncertain branches of a network. A subtype implements

- `uncertain_branches(model)`: the branch indices it overrides, and
- `branch_data(model, b, p, h, sequence)`: the per-unit π-model data
  `(br_r, br_x, g_fr, b_fr, g_to, b_to)` of branch `b` at the parameter point `p` (a named
  tuple), the harmonic order `h` (real) and the sequence `sequence` (`:positive`,
  `:negative` or `:zero`).

Branches it does not override keep the scaled data of [`replicate`](@ref).
"""
abstract type AbstractNetworkModel end

"""
    ScaledBranchModel(hdata, classes; lengths = Int[])

Analytic model for development and validation. Every uncertain branch `b` belongs to a
class `c = classes[b]` (a `Symbol`), and the parameters `kR_c`, `kL_c` and `kC_c` scale its
fundamental series resistance, series reactance and shunt susceptance. For the branches
in `lengths`, the parameter `len_b` scales all three, as a length factor. A parameter
absent from `p` is taken as 1. The harmonic scaling is that of [`replicate`](@ref):
r·√h, x·h, b·h. With every factor at 1, the data equal `replicate`'s exactly.

Positive and negative sequence only: the fundamental data carry no zero-sequence values.
"""
struct ScaledBranchModel <: AbstractNetworkModel
    classes::Dict{Int,Symbol}
    lengths::Set{Int}
    base::Dict{Int,NamedTuple}
end

function ScaledBranchModel(hdata::Dict, classes::Dict{Int,Symbol}; lengths = Int[])
    fund = hdata["nw"]["1"]["branch"]
    base = Dict{Int,NamedTuple}()
    for b in union(keys(classes), lengths)
        br = fund["$b"]
        base[b] = (br_r = br["br_r"], br_x = br["br_x"], g_fr = get(br, "g_fr", 0.0),
                   b_fr = get(br, "b_fr", 0.0), g_to = get(br, "g_to", 0.0),
                   b_to = get(br, "b_to", 0.0))
    end
    return ScaledBranchModel(classes, Set{Int}(lengths), base)
end

uncertain_branches(model::ScaledBranchModel) = sort(collect(keys(model.base)))

_factor(p, name::Symbol) = haskey(p, name) ? p[name] : 1.0

function branch_data(model::ScaledBranchModel, b::Int, p, h::Real, sequence::Symbol)
    sequence === :zero && throw(ArgumentError(
        "ScaledBranchModel: no zero-sequence data for branch $b; the fundamental data " *
        "are positive sequence"))
    c  = get(model.classes, b, nothing)
    kR = c === nothing ? 1.0 : _factor(p, Symbol("kR_", c))
    kL = c === nothing ? 1.0 : _factor(p, Symbol("kL_", c))
    kC = c === nothing ? 1.0 : _factor(p, Symbol("kC_", c))
    kℓ = b in model.lengths ? _factor(p, Symbol("len_", b)) : 1.0
    bd = model.base[b]
    return (br_r = kℓ * kR * bd.br_r * sqrt(h), br_x = kℓ * kL * bd.br_x * h,
            g_fr = kℓ * bd.g_fr, b_fr = kℓ * kC * bd.b_fr * h,
            g_to = kℓ * bd.g_to, b_to = kℓ * kC * bd.b_to * h)
end

"""
    ExactPiModel(branches; f1 = 50.0)

Model from per-unit-length data. `branches[b] = (zy, len, zbase)` gives, for branch `b`, a
function `zy(p, f, sequence) -> (z, y)` returning the series impedance `z` [Ω/m] and the
shunt admittance `y` [S/m] at frequency `f` [Hz], the length `len` [m] (a number, or a
function of `p` for an uncertain length), and the base impedance `zbase` [Ω]. The branch
data follow from the exact π-model of a uniform line: Z = Z_c·sinh(γℓ) and
Y/2 = tanh(γℓ/2)/Z_c, with γ = √(zy) and Z_c = √(z/y). This is the hook for cable and line
parameters computed by, e.g., LineCableModels.jl (see `examples/rhhc`).
"""
struct ExactPiModel <: AbstractNetworkModel
    branches::Dict{Int,NamedTuple}
    f1::Float64
end
ExactPiModel(branches::Dict{Int,<:NamedTuple}; f1 = 50.0) =
    ExactPiModel(Dict{Int,NamedTuple}(branches), f1)

uncertain_branches(model::ExactPiModel) = sort(collect(keys(model.branches)))

function branch_data(model::ExactPiModel, b::Int, p, h::Real, sequence::Symbol)
    br  = model.branches[b]
    len = br.len isa Function ? br.len(p) : br.len
    z, y = br.zy(p, model.f1 * h, sequence)
    γ, zc = sqrt(z * y), sqrt(z / y)
    zser  = zc * sinh(γ * len)
    yhalf = tanh(γ * len / 2) / zc
    return (br_r = real(zser) / br.zbase, br_x = imag(zser) / br.zbase,
            g_fr = real(yhalf) * br.zbase, b_fr = imag(yhalf) * br.zbase,
            g_to = real(yhalf) * br.zbase, b_to = imag(yhalf) * br.zbase)
end

# ──────────────────────────────────────────────────────────────────────────────
# Networks and admittance at a parameter point
# ──────────────────────────────────────────────────────────────────────────────

# Copy of `ntw` whose uncertain branches carry the data of `model` at (p, h, sequence).
# Only the branch table and the overridden branches are copied; everything else is shared.
function _override_branches(ntw::Dict, model::AbstractNetworkModel, p, h::Real,
                            sequence::Symbol)
    out = copy(ntw)
    out["branch"] = copy(ntw["branch"])
    for b in uncertain_branches(model)
        br = copy(ntw["branch"]["$b"])
        for (k, v) in pairs(branch_data(model, b, p, h, sequence))
            br[string(k)] = v
        end
        out["branch"]["$b"] = br
    end
    return out
end

"""
    networks_at(hdata, model, p) -> Dict

The harmonic networks of `hdata` (every nw id except the fundamental) with the uncertain
branches of `model` evaluated at the parameter point `p`, each at its own harmonic order and
sequence. Pass it to [`shhc_moments`](@ref) as `networks` to obtain the moments at `p`.
"""
function networks_at(hdata::Dict, model::AbstractNetworkModel, p)
    out = Dict{String,Any}()
    for (nw, ntw) in hdata["nw"]
        nw == "1" && continue
        h = parse(Int, nw)
        out[nw] = _override_branches(ntw, model, p, h, harmonic_sequence(h))
    end
    return out
end

"""
    harmonic_network(hdata, h; sequence = :positive) -> Dict

The network of `hdata` scaled to the (real) harmonic order `h` from its fundamental data,
with the rules of [`replicate`](@ref) (`_scale_network_to_harmonic!`). Between harmonics
the sequence is a modeling choice; positive sequence is the default.
"""
function harmonic_network(hdata::Dict, h::Real; sequence::Symbol = :positive)
    ntw = deepcopy(hdata["nw"]["1"])
    return _scale_network_to_harmonic!(ntw, h, sequence)
end

"""
    admittance_at(hdata, model, p, h; sequence = :positive) -> (Y, bus_idx)

Harmonic admittance matrix at the parameter point `p` and the (real) harmonic order `h`,
with the stamps of the stochastic second-order cone program (`_harmonic_admittance`), and
the map from bus id to matrix index. For the modal analysis that samples the box.
"""
function admittance_at(hdata::Dict, model::AbstractNetworkModel, p, h::Real;
                       sequence::Symbol = :positive)
    fund    = hdata["nw"]["1"]
    ids     = sort([bus["bus_i"] for (_, bus) in fund["bus"]])
    bus_idx = Dict(id => i for (i, id) in enumerate(ids))
    ntw     = _override_branches(harmonic_network(hdata, h; sequence), model, p, h, sequence)
    return _harmonic_admittance(ntw, bus_idx), bus_idx
end

# ──────────────────────────────────────────────────────────────────────────────
# Moments at a parameter point
# ──────────────────────────────────────────────────────────────────────────────

# The robust program needs a chance-constraint factor that does not depend on the network
# realization. Cantelli (distribution-free) and Gaussian are fixed numbers. The empirical
# factor λ* is calibrated per constraint on the sampled output of the NOMINAL network, so
# it does not carry over to another point of the box: reject it rather than reuse it.
function _check_robust_bound(hdata::Dict)
    pce   = hdata["pce"]
    bound = pce["bound"]
    bound in (:cantelli, :gaussian) || throw(ArgumentError(
        "robust HHC: the chance-constraint bound must be :cantelli or :gaussian, got " *
        ":$bound. The empirical factor λ* is calibrated on the nominal network and does " *
        "not hold at other points of the parameter box."))
    get(pce, "lambdas", nothing) === nothing || throw(ArgumentError(
        "robust HHC: per-constraint empirical factors λ* are present in hdata[\"pce\"]; " *
        "they are calibrated on the nominal network and cannot be used in the robust program."))
    return nothing
end

"""
    moments_at(hdata, model, p) -> NamedTuple

The moments of [`shhc_moments`](@ref) at the parameter point `p` of the network model
`model`. `hdata` must carry the PCE data with a Cantelli or Gaussian bound: the empirical
factor λ* is calibrated on the nominal network and is rejected.
"""
function moments_at(hdata::Dict, model::AbstractNetworkModel, p)
    _check_robust_bound(hdata)
    return shhc_moments(hdata; networks = networks_at(hdata, model, p))
end
