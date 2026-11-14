################################################################################
# HarmonicPowerModels.jl                                                       #
# Zero-sequence network of the second-order cone program (sHHC_SOC, RHHC).    #
#                                                                               #
# At a zero-sequence harmonic (h ≡ 0 mod 3) the network of the reduced         #
# program is not the positive-sequence one: transformers pass, earth or block  #
# zero-sequence current depending on their windings and neutral earthing,      #
# lines and cables have their own zero-sequence data, generators earth the     #
# network only through their zero-sequence impedance, and bus shunts are       #
# taken as unearthed stars. `_zero_sequence_view` converts a harmonic network  #
# into an equivalent network in the keys the positive-sequence code reads, so  #
# the admittance, the transfer impedances and the branch currents need no      #
# separate zero-sequence path. The IV formulations are not affected: they      #
# model the windings in their own constraints.                                  #
#                                                                               #
# Transformer windings (`vg` = "<side 1><side 2><clock>", e.g., "Yd11";        #
# `gnd1`/`gnd2` = 1 for an earthed star, neutral impedance re + j xe):          #
#   earthed–earthed       series   z + 3Zn1/|t|² + 3Zn2                        #
#   earthed–star point    series   z + 3Zn          (three-winding T model)     #
#   delta–earthed         shunt    z + 3Zn at the earthed side                  #
#   delta–star point      shunt    z at the star point                          #
#   any other pairing     open     (unearthed star, delta–delta, …)             #
# with z the leakage impedance (`br_r0 + j br_x0` if given, else the positive- #
# sequence `br_r + j br_x`). A star point is a bus with `star_point = true`.   #
################################################################################
# Author: Tom Van Acker                                                         #
################################################################################

import LinearAlgebra

# Configuration of side `s` (1 or 2) of a transformer winding pair: :star (the bus is the
# star point of a three-winding unit), :delta, :earthed or :open.
function _zs_side(cnf::Char, gnd, star::Bool)
    star && return :star
    cnf == 'D' && return :delta
    gnd == 1 && return :earthed
    return :open
end

# Zero-sequence equivalent of a winding pair: (:series, z), (:shunt_fr, z), (:shunt_to, z)
# or (:open, 0), for leakage z and neutral impedances zn1, zn2 (already ×3,
# with zn1 referred to the leakage side of any from-side tap).
function _zs_combination(s1::Symbol, s2::Symbol, z::Complex, zn1::Complex, zn2::Complex)
    s1 === :earthed && s2 === :earthed && return (:series, z + zn1 + zn2)
    s1 === :earthed && s2 === :star    && return (:series, z + zn1)
    s1 === :star    && s2 === :earthed && return (:series, z + zn2)
    s1 === :star    && s2 === :star    && return (:series, z)
    s1 === :earthed && s2 === :delta   && return (:shunt_fr, z + zn1)
    s1 === :delta   && s2 === :earthed && return (:shunt_to, z + zn2)
    s1 === :star    && s2 === :delta   && return (:shunt_fr, z)
    s1 === :delta   && s2 === :star    && return (:shunt_to, z)
    return (:open, zero(z))
end

# Zero-sequence polarity of two star windings, after phase relabelling.
_zs_star_angle(vg::AbstractString) = mod(parse(Int, vg[3:end]), 4) == 2 ? 180.0 : 0.0

_is_star_point(ntw::Dict, bus_id) =
    get(get(ntw["bus"], string(bus_id), Dict{String,Any}()), "star_point", false) == true

function _zs_transformer_branch!(br::Dict, ntw::Dict, id)
    vg = br["vg"]
    star1 = _is_star_point(ntw, br["f_bus"]); star2 = _is_star_point(ntw, br["t_bus"])
    for (w, star, cnf) in (("1", star1, uppercase(vg[1])), ("2", star2, uppercase(vg[2])))
        star || cnf == 'D' || haskey(br, "gnd" * w) || error(
            "zero sequence: transformer branch $id ($(vg)) has no gnd$w; declare the " *
            "earthing of its star winding (gnd$w = 1 earthed, 0 unearthed)")
    end
    s1 = _zs_side(uppercase(vg[1]), get(br, "gnd1", 0), star1)
    s2 = _zs_side(uppercase(vg[2]), get(br, "gnd2", 0), star2)
    z  = haskey(br, "br_x0") ? complex(get(br, "br_r0", 0.0), br["br_x0"]) :
                               complex(br["br_r"], br["br_x"])
    tm = get(br, "tap", 1.0); tm == 0.0 && (tm = 1.0)
    zn1 = 3 * complex(get(br, "re1", 0.0), get(br, "xe1", 0.0)) / abs2(tm)
    zn2 = 3 * complex(get(br, "re2", 0.0), get(br, "xe2", 0.0))
    kind, zeq = _zs_combination(s1, s2, z, zn1, zn2)
    br["br_r"], br["br_x"] = kind === :open ? (0.0, 0.0) : (real(zeq), imag(zeq))
    for k in ("g_fr", "b_fr", "g_to", "b_to"); br[k] = 0.0; end
    br["earth_fr"] = kind === :shunt_to          # the from side is earth (a delta winding)
    br["earth_to"] = kind === :shunt_fr
    # Star–star clocks 2, 6 and 10 reverse the polarity of all three windings.
    phase = kind === :series ? deg2rad(_zs_star_angle(vg)) : 0.0
    br["shift"]    = get(br, "shift", 0.0) - get(br, "vg_shift", 0.0) + phase
    br["vg_shift"] = phase
    return br
end

function _zs_plain_branch!(br::Dict, id)
    get(br, "sequence_data", nothing) === :zero && return br     # a network model's own data
    if get(br, "zero_sequence_open", false)
        br["br_r"], br["br_x"] = 0.0, 0.0
        for k in ("g_fr", "b_fr", "g_to", "b_to"); br[k] = 0.0; end
    elseif haskey(br, "br_x0")
        br["br_r"] = get(br, "br_r0", 0.0); br["br_x"] = br["br_x0"]
        br["g_fr"] = get(br, "g0_fr", 0.0); br["b_fr"] = get(br, "b0_fr", 0.0)
        br["g_to"] = get(br, "g0_to", 0.0); br["b_to"] = get(br, "b0_to", 0.0)
    else
        error("zero sequence: branch $id has no zero-sequence data; give br_r0, br_x0 " *
              "(and b0_fr, b0_to), model it with a network model that serves the zero " *
              "sequence, or set zero_sequence_open = true")
    end
    return br
end

"""
    _zero_sequence_view(ntw) -> Dict

The zero-sequence equivalent of the harmonic network `ntw` for the second-order cone
program, in the keys of the positive-sequence data (see the header of this file). Returns
`ntw` itself unless `ntw["harmonic_sequence"] === :zero`; the result is marked, so a second
call returns it unchanged. Missing zero-sequence data raise an error rather than fall back
to positive-sequence values.

- Branches with `vg`: the winding rules above, with `earth_fr`/`earth_to` for a shunt.
- Other branches: kept if a network model wrote zero-sequence data (`sequence_data`),
  else `br_r0`, `br_x0`, `g0_fr`, `b0_fr`, `g0_to`, `b0_to`, else open if
  `zero_sequence_open`, else an error.
- Generators: the zero-sequence admittance `gsc0 + j bsc0` (from `rsc0`, `xsc0`), none
  with `zero_sequence_open` or without a positive-sequence earth path, else an error.
- Bus shunts: unearthed stars (dropped) unless `earthed = true`.
- `xfmr` records: the winding rules, from `cnf1`/`cnf2`, `gnd1`/`gnd2`, `re`/`xe`, with
  a shunt added as an earthed bus shunt.
"""
function _zero_sequence_view(ntw::Dict)
    get(ntw, "harmonic_sequence", :positive) === :zero || return ntw
    get(ntw, "zero_sequence_view", false) && return ntw
    out = copy(ntw)
    out["zero_sequence_view"] = true

    out["branch"] = Dict{String,Any}()
    for (id, br0) in ntw["branch"]
        br = copy(br0)
        haskey(br, "vg") ? _zs_transformer_branch!(br, ntw, id) : _zs_plain_branch!(br, id)
        out["branch"][id] = br
    end

    out["gen"] = Dict{String,Any}()
    for (id, g0) in ntw["gen"]
        g = copy(g0)
        if !iszero(get(g, "inf", 0))
            # an ideal source: no shunt in either sequence
        elseif haskey(g, "gsc0")
            g["gsc"], g["bsc"] = g["gsc0"], g["bsc0"]
        elseif get(g, "zero_sequence_open", false) ||
               (iszero(get(g, "gsc", 0.0)) && iszero(get(g, "bsc", 0.0)))
            g["gsc"], g["bsc"] = 0.0, 0.0
        else
            error("zero sequence: generator $id has no zero-sequence impedance; give rsc0, " *
                  "xsc0 (e.g., the leakage of its YNd step-up) or set zero_sequence_open = true")
        end
        out["gen"][id] = g
    end

    out["shunt"] = Dict{String,Any}(id => s for (id, s) in get(ntw, "shunt", Dict{String,Any}())
                                    if get(s, "earthed", false))

    if haskey(ntw, "xfmr")
        out["xfmr"] = Dict{String,Any}()
        for (id, x0) in ntw["xfmr"]
            x  = copy(x0)
            s1 = _zs_side(Char(get(x, "cnf1", 'Y')), get(x, "gnd1", 0), false)
            s2 = _zs_side(Char(get(x, "cnf2", 'Y')), get(x, "gnd2", 0), false)
            z  = complex(get(x, "r1", 0.0) + get(x, "r2", 0.0), get(x, "xsc", 0.0))
            kind, zeq = _zs_combination(s1, s2, z,
                3 * complex(get(x, "re1", 0.0), get(x, "xe1", 0.0)),
                3 * complex(get(x, "re2", 0.0), get(x, "xe2", 0.0)))
            if kind === :series
                x["r1"], x["r2"], x["xsc"] = real(zeq), 0.0, imag(zeq)
                angle = _zs_star_angle(get(x, "vg", "Yy0"))
                x["tr"], x["ti"] = cosd(angle), sind(angle)
                out["xfmr"][id] = x
            elseif kind !== :open
                bus = kind === :shunt_fr ? x["f_bus"] : x["t_bus"]
                y   = 1 / zeq
                out["shunt"]["xfmr_zs_$id"] = Dict{String,Any}("shunt_bus" => bus, "gs" => real(y),
                    "bs" => imag(y), "earthed" => true, "status" => 1)
            end
        end
    end
    return out
end

"""
    _earth_floating(Y; tol = 1e-12) -> (Y′, floating)

For a reduced zero-sequence admittance `Y`, the connected parts with no path to earth
(singular blocks, e.g., behind delta windings without any shunt) and a copy of `Y` with the
identity on their diagonal, which makes it invertible without changing the other parts.
Units in a floating part are omitted from the injections and allocation weights at
this harmonic (see `shhc_moments`); units in earthed parts retain their capacity. Only
when every harmonic unit floats is the harmonic blocked with zero capacity. Voltages
in floating parts are zero for injections elsewhere.
Returns `Y` itself and an empty list when no part floats.
"""
function _earth_floating(Y::AbstractMatrix; tol::Real = 1e-12)
    n = size(Y, 1)
    comp = zeros(Int, n); nc = 0
    for s in 1:n
        comp[s] == 0 || continue
        nc += 1; comp[s] = nc; stack = [s]
        while !isempty(stack)
            u = pop!(stack)
            for v in 1:n
                (comp[v] == 0 && (Y[u, v] != 0 || Y[v, u] != 0)) || continue
                comp[v] = nc; push!(stack, v)
            end
        end
    end
    floating = Vector{Int}[]
    for c in 1:nc
        idx = findall(==(c), comp)
        B   = Y[idx, idx]
        sv  = LinearAlgebra.svdvals(Matrix(B))
        (isempty(sv) || maximum(sv) == 0 || minimum(sv) <= tol * maximum(sv)) && push!(floating, idx)
    end
    isempty(floating) && return Y, floating
    Y2 = Matrix(copy(Y))
    for idx in floating, i in idx
        Y2[i, i] += 1.0
    end
    return Y2, floating
end
