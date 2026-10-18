################################################################################
# HarmonicPowerModels.jl                                                       #
# Harmonic admittance matrix and transfer-impedance utilities.                 #
#                                                                               #
# Three public functions:                                                       #
#   build_harmonic_Ybus             – full n_bus × n_bus nodal admittance      #
#   harmonic_transfer_impedance     – Z mapping load currents → bus voltages   #
#   bus_voltage_limits              – squared IHD voltage limits per bus       #
################################################################################
# Authors: Tom Van Acker                                                        #
################################################################################

using LinearAlgebra: factorize, I

# ──────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ──────────────────────────────────────────────────────────────────────────────

function _bus_id_map(data::Dict)
    ids = sort([parse(Int, k) for k in keys(data["bus"])])
    return ids, Dict(id => i for (i, id) in enumerate(ids))
end

# ──────────────────────────────────────────────────────────────────────────────
# 1. build_harmonic_Ybus
# ──────────────────────────────────────────────────────────────────────────────

"""
    build_harmonic_Ybus(data, h) -> Y::Matrix{ComplexF64}

Assemble the (n_bus × n_bus) nodal admittance matrix for harmonic order `h`.

Bus IDs are sorted from `data["bus"]` and mapped to 1-indexed Julia rows/cols.
The fundamental-frequency `data` values are scaled to harmonic `h` internally:

  - branch:  br_r·√h,  br_x·h,  b_fr·h,  b_to·h
  - shunt:   gs (frequency-independent),  bs·h if capacitive (bs ≥ 0),
             bs/h if inductive (bs < 0) — matching `replicate`
  - xfmr:    xsc·h  (r1, r2 ignored — lossless model)
  - xfmr phase shift: pos-seq (h%3==1) ∠−30·shift°,
                      neg-seq (h%3==2) ∠+30·shift°,
                      zero-seq (h%3==0) ∠0°

**Zero-sequence delta blocking** (h%3==0): for each transformer whose `vg`
string has 'D' (delta) on the primary or secondary side, the corresponding bus
row and column are zeroed after all stamps, modelling the isolation of that bus
from zero-sequence currents.

The magnetising branch (gsh) and winding resistances (r1, r2) are not included.
"""
function build_harmonic_Ybus(data::Dict, h::Int)
    bus_ids, bus_idx = _bus_id_map(data)
    n_bus = length(bus_ids)

    Y = zeros(ComplexF64, n_bus, n_bus)

    # ── branches (π-model) ────────────────────────────────────────────────────
    for branch in values(data["branch"])
        get(branch, "br_status", 1) == 0 && continue

        fr = bus_idx[branch["f_bus"]]
        to = bus_idx[branch["t_bus"]]

        r_h    = branch["br_r"]              * sqrt(h)
        x_h    = branch["br_x"]              * h
        b_fr_h = get(branch, "b_fr", 0.0)   * h   # b_fr = b_c/2 in PowerModels
        b_to_h = get(branch, "b_to", 0.0)   * h

        y_s  = 1 / (r_h + im * x_h)
        y_fr = im * b_fr_h
        y_to = im * b_to_h

        Y[fr, fr] += y_s + y_fr
        Y[fr, to] -= y_s
        Y[to, fr] -= y_s
        Y[to, to] += y_s + y_to
    end

    # ── bus shunts (capacitor banks, shunt filters, shunt reactors) ───────────
    # Same frequency scaling as `replicate`: gs is frequency-independent, a
    # capacitive bs (≥ 0) scales with h, an inductive bs (< 0) with 1/h.
    for shunt in values(get(data, "shunt", Dict{String,Any}()))
        get(shunt, "status", 1) == 0 && continue
        s    = bus_idx[shunt["shunt_bus"]]
        bs   = get(shunt, "bs", 0.0)
        bs_h = bs >= 0.0 ? bs * h : bs / h
        Y[s, s] += get(shunt, "gs", 0.0) + im * bs_h
    end

    # ── transformers ──────────────────────────────────────────────────────────
    if haskey(data, "xfmr")
        is_zero_seq  = (h % 3 == 0)
        blocked      = Set{Int}()

        for xfmr in values(data["xfmr"])
            fr = bus_idx[xfmr["f_bus"]]
            to = bus_idx[xfmr["t_bus"]]

            y_sc = 1 / (im * xfmr["xsc"] * h)

            vg    = get(xfmr, "vg", "Yy0")
            cnf1  = uppercase(vg[1])
            cnf2  = uppercase(vg[2])
            shift = parse(Int, vg[3:end])

            if h % 3 == 1        # positive-sequence
                tr = cosd(-30.0 * shift);  ti = sind(-30.0 * shift)
            elseif h % 3 == 2    # negative-sequence
                tr = cosd( 30.0 * shift);  ti = sind( 30.0 * shift)
            else                 # zero-sequence
                tr = 1.0;                  ti = 0.0
            end
            t = tr + im * ti

            # Stamp (valid for |t| = 1 phase-shifting transformers)
            Y[fr, fr] += y_sc
            Y[fr, to] -= y_sc * t
            Y[to, fr] -= conj(t) * y_sc
            Y[to, to] += y_sc

            # Collect delta-blocked buses for zero-sequence
            if is_zero_seq
                cnf1 == 'D' && push!(blocked, fr)
                cnf2 == 'D' && push!(blocked, to)
            end
        end

        # Zero rows and columns of delta-blocked buses
        for b in blocked
            Y[b, :] .= 0
            Y[:, b] .= 0
        end
    end

    return Y
end

# ──────────────────────────────────────────────────────────────────────────────
# 2. harmonic_transfer_impedance
# ──────────────────────────────────────────────────────────────────────────────

"""
    harmonic_transfer_impedance(data, h, load_bus_ids) -> Z::Matrix{ComplexF64}

Compute the transfer-impedance sub-matrix for harmonic `h`:

    Z[:, u]  =  V_active   when bus load_bus_ids[u] injects 1+0j A
                            and all other buses inject 0 A

Returns a matrix of shape `(n_active × n_loads)` where `n_active` is the
number of buses remaining after Kron-reduction of the slack bus (bus_type=3)
and any zero-sequence delta-blocked buses (buses with all-zero Y rows).  The
row ordering follows the sorted active bus IDs.

For positive/negative-sequence harmonics (h%3 ≠ 0) there is no delta blocking
and `n_active = n_bus − 1`.
"""
function harmonic_transfer_impedance(data::Dict, h::Int,
                                      load_bus_ids::Vector{Int})
    bus_ids, bus_idx = _bus_id_map(data)
    n_bus = length(bus_ids)

    # Locate slack bus (bus_type == 3)
    slack_mat = 0
    for (k, bus) in data["bus"]
        if get(bus, "bus_type", 1) == 3
            slack_mat = bus_idx[parse(Int, k)]
            break
        end
    end
    slack_mat > 0 || error("build_harmonic_Ybus: no slack bus (bus_type=3) in data")

    Y_full = build_harmonic_Ybus(data, h)

    # Buses to eliminate: slack + any bus whose Y row is entirely zero
    # (delta-blocked buses in zero-sequence, isolated buses, etc.)
    elim = Set{Int}([slack_mat])
    for i in 1:n_bus
        all(iszero, Y_full[i, :]) && push!(elim, i)
    end

    active  = [i for i in 1:n_bus if i ∉ elim]
    n_act   = length(active)
    n_loads = length(load_bus_ids)
    Z       = zeros(ComplexF64, n_act, n_loads)

    n_act == 0 && return Z

    Y_act   = Y_full[active, active]
    act_idx = Dict(mat => pos for (pos, mat) in enumerate(active))

    Yfact = factorize(Y_act)

    for (u, lb_id) in enumerate(load_bus_ids)
        mat = bus_idx[lb_id]
        haskey(act_idx, mat) || continue   # load bus eliminated (e.g. delta-blocked)
        I_unit             = zeros(ComplexF64, n_act)
        I_unit[act_idx[mat]] = 1.0 + 0.0im
        Z[:, u] = Yfact \ I_unit
    end

    return Z
end

# ──────────────────────────────────────────────────────────────────────────────
# 3. bus_voltage_limits
# ──────────────────────────────────────────────────────────────────────────────

"""
    bus_voltage_limits(data, h) -> limits::Dict{Int, Float64}

Return the squared IHD voltage limit at each bus for harmonic `h`:

    limits[bus_id]  =  (ihdmax_{bus,h} · Vfund_{bus})²

where `ihdmax` comes from `ihd_limits[bus["standard"]][h]` and `Vfund` is
the fundamental voltage magnitude `bus["vm"]`.

!!! note "Fundamental voltage reference"
    `Vfund = bus["vm"]` is taken as-is from `data`. On raw parsed data this is
    the **nominal** setpoint (typically 1.0), so the IHD limit is relative to
    nominal voltage. The JuMP solvers (`sHHC_QCQP`, `sHHC_SOC`) instead use the
    **HPF-solved** fundamental voltage, so on networks with appreciable voltage
    drop the solver limit is tighter than this one. To make this function
    consistent with the solvers, populate `bus["vm"]` with HPF results before
    calling (e.g. via `update_hdata_with_fundamental_hpf_results!`).

Buses without a `standard` key, or whose standard does not include harmonic
`h`, receive a limit of 0.0 (no allowance).
"""
function bus_voltage_limits(data::Dict, h::Int)
    limits = Dict{Int, Float64}()
    for (k, bus) in data["bus"]
        bus_id = parse(Int, k)
        vm     = get(bus, "vm", 1.0)

        ihdmax = 0.0
        if haskey(bus, "standard")
            std = bus["standard"]
            if haskey(ihd_limits, std) && 1 <= h <= length(ihd_limits[std])
                ihdmax = ihd_limits[std][h]
            end
        end

        limits[bus_id] = (ihdmax * vm)^2
    end
    return limits
end

# ──────────────────────────────────────────────────────────────────────────────
# 4. Constructive worst-case angle alignment (SHHC paper baseline)
# ──────────────────────────────────────────────────────────────────────────────

"""
    reduced_harmonic_impedance(hdata) -> (Z, bus_pos)

Kron-reduced harmonic transfer impedance for every harmonic network in the
multinetwork `hdata` (nw ≠ "1").

For each harmonic order `h`, the nodal admittance is built with
`_harmonic_admittance` from the **already harmonic-scaled** network data, the
reference (slack, `bus_type == 3`) buses are eliminated, and the reduced matrix
is inverted to obtain the transfer impedance between non-reference buses. This
is the same `Z_h` used internally by `solve_shhc(::sHHC_SOC)`, so the resulting
angles are exactly consistent with the SOCP model.

Returns
  - `Z::Dict{Int,Matrix{ComplexF64}}` keyed by harmonic order; row/column `k`
    corresponds to the non-reference bus with id `nonref_ids[k]`.
  - `bus_pos::Dict{Int,Int}` mapping bus id → matrix index (non-reference buses).

`Z[h][bus_pos[n], bus_pos[d]]` is the harmonic voltage at bus `n` per unit
current injected at bus `d`.
"""
function reduced_harmonic_impedance(hdata::Dict)
    fund    = hdata["nw"]["1"]
    bus_ids = sort([b["bus_i"] for (_, b) in fund["bus"]])
    bus_idx = Dict(id => i for (i, id) in enumerate(bus_ids)); Nn = length(bus_ids)
    # Clean buses (bus_type 4) are ideal harmonic grounds and are Kron-eliminated;
    # every other bus is grounded through the generator subtransient admittance
    # Y_g(h) (stamped in _harmonic_admittance) and the line/bus shunts, keeping
    # the reduced matrix invertible (consistent with solve_shhc(::sHHC_SOC)).
    refset  = Set(bus_idx[b["bus_i"]] for (_, b) in fund["bus"] if get(b, "is_clean", false))
    nonref  = [i for i in 1:Nn if i ∉ refset]
    nonref_ids = [bus_ids[i] for i in nonref]
    bus_pos    = Dict(id => k for (k, id) in enumerate(nonref_ids))

    Z = Dict{Int, Matrix{ComplexF64}}()
    for (nw_str, ntw) in hdata["nw"]
        nw_str == "1" && continue
        h  = parse(Int, nw_str)
        Yh = _harmonic_admittance(ntw, bus_idx)
        Yr = Yh[nonref, nonref]
        try
            Z[h] = Yr \ Matrix{ComplexF64}(I, length(nonref), length(nonref))
        catch
            @warn "reduced_harmonic_impedance: singular reduced admittance for " *
                  "harmonic $h; skipped."
        end
    end
    return Z, bus_pos
end

"""
    constructive_worst_case_angles(Z, bus_pos, hdata, nstar) -> Dict{Int,Dict{Int,Float64}}

Per-harmonic, per-unit harmonic-current angle (**radians**) that aligns every
unit's contribution IN PHASE at the binding bus `nstar`:

    φ_{h,d} = −arg( Z_h[nstar, bus(d)] )

Because each harmonic is aligned independently, this maximises Σ_h |V_{nstar,h}|²
— the THD / RMS-voltage at `nstar` — and therefore yields the maximum harmonic
distortion (the LOWER BOUND of hosting capacity). `Z`/`bus_pos` come from
`reduced_harmonic_impedance`. Every load gets an entry; units whose bus is
eliminated (reference / delta-blocked — no harmonic-voltage coupling, so the
angle is irrelevant) receive `0.0`. Returns `Dict[h][load_id] = angle` (radians).
"""
function constructive_worst_case_angles(Z::Dict, bus_pos::Dict, hdata::Dict,
                                        nstar::Int)
    haskey(bus_pos, nstar) || error("constructive_worst_case_angles: bus " *
                                    "$nstar is not a non-reference bus")
    np = bus_pos[nstar]
    angles = Dict{Int, Dict{Int, Float64}}()
    for (nw_str, ntw) in hdata["nw"]
        nw_str == "1" && continue
        h = parse(Int, nw_str)
        haskey(Z, h) || continue
        ah = Dict{Int, Float64}()
        for (d_str, load) in ntw["load"]
            dp = get(bus_pos, load["load_bus"], 0)
            ah[parse(Int, d_str)] = dp == 0 ? 0.0 : -angle(Z[h][np, dp])
        end
        angles[h] = ah
    end
    return angles
end

"""
    empirical_lambdas(hdata, g; nsamp=20000, seed=1) -> Dict{Symbol,Dict}

Per-constraint empirical (quantile-calibrated) chance-constraint factors

    λ*(ε) = ( quantile_{1−ε}(S) − E[S] ) / σ(S)

estimated by Monte-Carlo sampling the OUTPUT distribution of each constrained
aggregate `S`, using the exact reduced harmonic transfer impedance and the
per-unit angle supports stored in `hdata["pce"]`. One Beta germ per unit
(independent across units, comonotone across that unit's harmonics).

`g` maps harmonic order → I_h² (the current hosting decision); it enters only
the cross-harmonic aggregates (THD / RMS), so the IHD factors are g-independent.

Returns `Dict(:ihd => Dict((h,bus)=>λ), :agg => Dict(bus=>λ), :rmsi => Dict(branch=>λ))`.
For `S = a·g_h` (per-harmonic IHD) `λ*` does not depend on `g`; for the
cross-harmonic aggregates it is evaluated at the supplied `g` (⇒ the outer
solve→recalibrate iteration). Used with `pce_params["bound"]="empirical"` and
`pce_params["empirical_lambdas"]=…`.
"""
function empirical_lambdas(hdata::Dict, g::AbstractDict; nsamp::Int=20000, seed::Int=1)
    pce    = hdata["pce"]
    ε      = pce["epsilon"]
    α̂, β̂   = pce["alpha_beta"]
    ranges = pce["angle_ranges"]
    units  = pce["harmonic_units"]
    Z, bus_pos = reduced_harmonic_impedance(hdata)
    fund   = hdata["nw"]["1"]
    ubus   = Dict(d => fund["load"]["$d"]["load_bus"] for d in units if haskey(fund["load"], "$d"))
    harm_nws = sort([parse(Int, s) for s in keys(hdata["nw"]) if s != "1" && haskey(Z, parse(Int, s))])

    npos    = length(bus_pos)
    pos2bus = Dict(p => n for (n, p) in bus_pos)

    # branch endpoints + series admittance per harmonic (positions; 0 = reference)
    brinfo = Dict(h => Tuple{Int,Int,Int,Float64,Float64}[] for h in harm_nws)
    for h in harm_nws
        for (b_str, br) in hdata["nw"]["$h"]["branch"]
            r, x = br["br_r"], br["br_x"]; (abs(r) + abs(x) > 0.0) || continue
            ys = 1.0 / (r + im * x)
            fp = get(bus_pos, br["f_bus"], 0); tp = get(bus_pos, br["t_bus"], 0)
            (fp == 0 && tp == 0) && continue
            push!(brinfo[h], (parse(Int, b_str), fp, tp, real(ys), imag(ys)))
        end
    end

    aV = Dict(h => [Vector{Float64}(undef, nsamp) for _ in 1:npos] for h in harm_nws)
    aI = Dict(h => Dict(t[1] => Vector{Float64}(undef, nsamp) for t in brinfo[h]) for h in harm_nws)

    rng = Random.MersenneTwister(seed); Bd = Distributions.Beta(α̂, β̂)
    Vr = zeros(npos); Vi = zeros(npos)
    ωu = Dict(d => 0.0 for d in units)
    for m in 1:nsamp
        for d in units; ωu[d] = rand(rng, Bd); end
        for h in harm_nws
            fill!(Vr, 0.0); fill!(Vi, 0.0); Zh = Z[h]
            for d in units
                haskey(ubus, d) && haskey(bus_pos, ubus[d]) || continue
                dp = bus_pos[ubus[d]]; lo, hi = ranges[h][d]
                φ = lo + (hi - lo) * ωu[d]; c = cos(φ); s = sin(φ)
                @inbounds for p in 1:npos
                    z = Zh[p, dp]; zr = real(z); zi = imag(z)
                    Vr[p] += zr * c - zi * s; Vi[p] += zi * c + zr * s
                end
            end
            @inbounds for p in 1:npos; aV[h][p][m] = Vr[p]^2 + Vi[p]^2; end
            for (bid, fp, tp, gy, by) in brinfo[h]
                dVr = (fp == 0 ? 0.0 : Vr[fp]) - (tp == 0 ? 0.0 : Vr[tp])
                dVi = (fp == 0 ? 0.0 : Vi[fp]) - (tp == 0 ? 0.0 : Vi[tp])
                Ir = gy * dVr - by * dVi; Ii = gy * dVi + by * dVr
                aI[h][bid][m] = Ir^2 + Ii^2
            end
        end
    end

    λstar(x) = (σ = Statistics.std(x); σ > 0 ? max((Statistics.quantile(x, 1 - ε) -
                Statistics.mean(x)) / σ, 0.0) : 0.0)

    λ_ihd = Dict{Tuple{Int,Int},Float64}()
    for h in harm_nws, p in 1:npos
        λ_ihd[(h, pos2bus[p])] = λstar(aV[h][p])
    end
    gv(h) = get(g, h, 0.0)
    λ_agg = Dict{Int,Float64}()
    for p in 1:npos
        S = sum(aV[h][p] .* gv(h) for h in harm_nws)
        λ_agg[pos2bus[p]] = λstar(S)
    end
    λ_rmsi = Dict{Int,Float64}()
    allbr = union((Set(keys(aI[h])) for h in harm_nws)...)
    for bid in allbr
        S = sum((haskey(aI[h], bid) ? aI[h][bid] : zeros(nsamp)) .* gv(h) for h in harm_nws)
        λ_rmsi[bid] = λstar(S)
    end
    return Dict(:ihd => λ_ihd, :agg => λ_agg, :rmsi => λ_rmsi)
end

"""
    recentre_angle_ranges(base, centres, α̂, β̂) -> Dict{Int,Dict{Int,Tuple}}

Re-centre each harmonic's Beta angle spread on a per-unit centre angle so that
its expectation equals that centre.

`base[h] = (v_min, v_max)` is the section-II angle range for harmonic `h`
(radians); `centres[h][d]` is the per-unit centre (radians, e.g. the constructive
worst-case angle). With Beta mean fraction `Eω = α̂/(α̂+β̂)`, the section-II mean is
`Eᵇ = v_min + (v_max−v_min)·Eω`, and the recentred per-unit range is

    (lo, hi) = ( centre + (v_min − Eᵇ),  centre + (v_max − Eᵇ) )

so that `E[φ_{h,d}] = centre` while the spread/shape are preserved. The output is
suitable for `pce_params["harmonics"]` (per-unit `Dict[h][d] = (lo,hi)`).
"""
function recentre_angle_ranges(base::Dict, centres::Dict{Int,<:AbstractDict},
                               α̂::Real, β̂::Real)
    Eω  = α̂ / (α̂ + β̂)
    out = Dict{Int, Dict{Int, Tuple{Float64,Float64}}}()
    for (h, c_units) in centres
        haskey(base, h) || continue
        vmin, vmax = base[h]
        out[h] = Dict(d => recentre_range(vmin, vmax, c, Eω) for (d, c) in c_units)
    end
    return out
end
