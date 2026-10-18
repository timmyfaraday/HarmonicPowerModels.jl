################################################################################
# HarmonicPowerModels.jl                                                       #
# LP reduction of the absolute-equality stochastic HHC problem.               #
#                                                                               #
# Under absolute-equality fairness all loads at harmonic h receive the same   #
# scalar hosting capacity I_h.  Bus voltages are then linear in I_h, and the  #
# Cantelli IHD chance constraint reduces to the scalar bound                  #
#                                                                               #
#   I_h ≤  sqrt( V_lim²_{n,h} / C_{n,h} )                                    #
#                                                                               #
# where, in the probability-measure PCE basis (see core/pce_basis.jl),        #
#   C_{n,h} = alpha_{n,h,0} + λ sqrt(Σ_{k≠0} ⟨Ψ_k,Ψ_k⟩ · alpha_{n,h,k}²) #
# and alpha_{n,h,k} = Σ_{κ1,κ2} M[κ1,κ2,k] (C_re[n,κ1]C_re[n,κ2]           #
#                                            + C_im[n,κ1]C_im[n,κ2])           #
# with C_re, C_im being the real and imaginary voltage-response PCE matrices.  #
#                                                                               #
# Reference: Section IV.D of the SHHC paper.                                  #
# Authors: Tom Van Acker                                                        #
################################################################################

################################################################################
# Internal helpers                                                              #
################################################################################

"""
    _active_bus_ids_lp(data, h) -> Vector{Int}

Sorted bus IDs that remain after Kron-elimination of the slack bus and any
zero-sequence delta-blocked buses (rows with all-zero Y entries) for harmonic h.
These correspond 1-to-1 with the rows of the Z_h matrix returned by
`harmonic_transfer_impedance`.
"""
function _active_bus_ids_lp(data::Dict, h::Int)
    bus_ids, bus_idx = _bus_id_map(data)
    n_bus = length(bus_ids)

    slack_mat = 0
    for (k, bus) in data["bus"]
        if get(bus, "bus_type", 1) == 3
            slack_mat = bus_idx[parse(Int, k)]
            break
        end
    end
    slack_mat > 0 || error("_active_bus_ids_lp: no slack bus found")

    Y_full = build_harmonic_Ybus(data, h)

    elim = Set{Int}([slack_mat])
    for i in 1:n_bus
        all(iszero, Y_full[i, :]) && push!(elim, i)
    end

    return bus_ids[[i for i in 1:n_bus if i ∉ elim]]
end

################################################################################
# 1. compute_C_coefficients                                                     #
################################################################################

"""
    compute_C_coefficients(Z_h, cos_coefs, sin_coefs, harmonics_units, pb, h)
        -> (C_re, C_im)::Tuple{Matrix{Float64}, Matrix{Float64}}
                            # shape (n_active_buses × pb.dim)

Real and imaginary voltage-response PCE coefficient matrices for harmonic `h`.

For each active bus `n` (row of Z_h) and PCE index k ∈ 1:pb.dim:

    C_re[n,k] = Σ_u [ Re(Z_h[n,u]) · cos_k[u]  −  Im(Z_h[n,u]) · sin_k[u] ]
    C_im[n,k] = Σ_u [ Im(Z_h[n,u]) · cos_k[u]  +  Re(Z_h[n,u]) · sin_k[u] ]

where `u` runs over the columns of Z_h (one per load unit), and:

    cos_k[u]  =  cos_coefs[(harmonics_units[u], h)][k]
    sin_k[u]  =  sin_coefs[(harmonics_units[u], h)][k]

(1-indexed Julia arrays; k=1 is the constant / zero-degree PCE coefficient.)

# Arguments
- `Z_h`            — (n_active × n_loads) complex transfer-impedance matrix.
- `cos_coefs`      — Dict{Tuple{Any,Int},Vector{Float64}} from `angle_pce_coefficients`.
- `sin_coefs`      — same.
- `harmonics_units`— load unit IDs, one per column of Z_h (same order).
- `pb`             — `PCEBasis` the coefficient vectors live in.
- `h`              — harmonic order for dict key lookup.
"""
function compute_C_coefficients(Z_h          ::Matrix{ComplexF64},
                                 cos_coefs   ::Dict,
                                 sin_coefs   ::Dict,
                                 harmonics_units::AbstractVector,
                                 pb          ::PCEBasis,
                                 h           ::Int)
    n_buses = size(Z_h, 1)
    P       = pb.dim

    C_re = zeros(n_buses, P)
    C_im = zeros(n_buses, P)

    for (u_col, u_id) in enumerate(harmonics_units)
        cos_k = cos_coefs[(u_id, h)]   # length-P vector
        sin_k = sin_coefs[(u_id, h)]
        for n in 1:n_buses
            rz = real(Z_h[n, u_col])
            iz = imag(Z_h[n, u_col])
            for k in 1:P
                C_re[n, k] += rz * cos_k[k] - iz * sin_k[k]
                C_im[n, k] += iz * cos_k[k] + rz * sin_k[k]
            end
        end
    end

    return C_re, C_im
end

################################################################################
# 2. compute_alpha                                                              #
################################################################################

"""
    compute_alpha(C_re, C_im, pb::PCEBasis) -> alpha::Matrix{Float64}
                                                 # shape (n_buses × pb.dim)

PCE coefficients of the squared voltage magnitude |V|²/I², per bus, via the
sparse Galerkin product:

    alpha[n, :] = C_re[n,:] ⊛ C_re[n,:]  +  C_im[n,:] ⊛ C_im[n,:]

where ⊛ is `galerkin_product`. Because the basis is orthogonal w.r.t. the
probability measure, alpha[n, 1] = E[|V|²/I²] directly.
"""
function compute_alpha(C_re::Matrix{Float64},
                       C_im::Matrix{Float64},
                       pb  ::PCEBasis)
    n_buses = size(C_re, 1)
    alpha   = zeros(n_buses, pb.dim)

    for n in 1:n_buses
        alpha[n, :] .= galerkin_product(pb, view(C_re, n, :), view(C_re, n, :)) .+
                       galerkin_product(pb, view(C_im, n, :), view(C_im, n, :))
    end

    # alpha[:,1] = E[|V|²/I²] must be non-negative
    @assert all(alpha[n, 1] >= -1e-10 for n in 1:n_buses) "alpha[n, k=0] must be non-negative for all buses"

    return alpha
end

################################################################################
# 3. compute_C_scalar                                                           #
################################################################################

"""
    compute_C_scalar(alpha, pb::PCEBasis, epsilon) -> C_scalar::Vector{Float64}
                                                       # length n_buses

Cantelli coefficient per bus from the PCE of |V|²/I² (probability-measure
basis, so no total-measure correction is needed):

    lambda      = √((1 − ε) / ε)
    E[|V|²/I²]  = alpha[n, 1]
    Var[|V|²/I²]= Σ_{k=2}^{P} ⟨Ψ_k,Ψ_k⟩ · alpha[n,k]²
    C_scalar[n] = E + lambda · √Var

C_scalar[n] is the denominator of the IHD chance constraint:
    I_h ≤ sqrt( V_lim²[n] / C_scalar[n] )
"""
function compute_C_scalar(alpha  ::Matrix{Float64},
                           pb     ::PCEBasis,
                           epsilon::Float64;
                           bound  ::Symbol = :cantelli)
    n_buses  = size(alpha, 1)
    lambda   = risk_factor(epsilon, bound)
    C_scalar = zeros(n_buses)

    for n in 1:n_buses
        var_n       = sum(pb.sp2[k] * alpha[n, k]^2 for k in 2:pb.dim; init=0.0)
        C_scalar[n] = alpha[n, 1] + lambda * sqrt(var_n)
    end

    @assert all(C_scalar[n] > -1e-10 for n in 1:n_buses) "C_scalar must be positive for all active buses"

    return C_scalar
end

################################################################################
# 4. lp_solution                                                                #
################################################################################

"""
    lp_solution(C_scalar, voltage_limits_sq, bus_indices)
        -> (I_star::Float64, binding_bus::Int)

Solve the single-harmonic LP reduction:

    I_star = min_{n : V_lim²[n] > 0}  sqrt( V_lim²[n] / C_scalar[n] )

where the minimum is taken only over buses that have a positive voltage limit.
Returns the optimal hosting capacity `I_star` and the ID of the binding bus.

# Arguments
- `C_scalar`           — length-n_active vector from `compute_C_scalar`.
- `voltage_limits_sq`  — Dict{Int,Float64} from `bus_voltage_limits`; keyed by bus ID.
- `bus_indices`        — bus IDs corresponding to rows of `C_scalar` (same order).
"""
function lp_solution(C_scalar         ::AbstractVector{Float64},
                     voltage_limits_sq ::Dict{Int, Float64},
                     bus_indices       ::AbstractVector{Int})
    I_star      = Inf
    binding_bus = -1

    for (n, bus_id) in enumerate(bus_indices)
        lim_sq = get(voltage_limits_sq, bus_id, 0.0)
        lim_sq > 0.0       || continue   # no harmonic limit at this bus
        C_scalar[n] > 0.0  || continue   # degenerate: bus unreachable from loads

        I_cand = sqrt(lim_sq / C_scalar[n])
        if I_cand < I_star
            I_star      = I_cand
            binding_bus = bus_id
        end
    end

    return I_star, binding_bus
end

################################################################################
# 5. compare_shared_vs_independent                                              #
################################################################################

"""
    compare_shared_vs_independent(data, h, epsilon, beta_params, angle_params)
        -> NamedTuple

Run the LP reduction twice and compare the optimal hosting capacity under
two correlation models:

(a) **Independent germs** — each load unit has its own Beta-distributed germ;
    germ space has dimension `n_loads`.
(b) **Shared germ** — all load units are driven by the same Beta germ;
    germ space has dimension 1.  All units receive the same cos/sin PCE
    coefficient vectors (from the first unit's angle bounds), and C_re[n,k]
    is the sum over all units for each k.

# Arguments
- `data`          — PowerModels-parsed network data dict.
- `h`             — harmonic order.
- `epsilon`       — Cantelli risk level ε ∈ (0,1).
- `beta_params`   — `NamedTuple` with fields:
    - `alpha_shape::Float64` — Beta shape parameter α̂
    - `beta_shape::Float64`  — Beta shape parameter β̂
    - `delta::Int`           — PCE degree
    - `n_terms::Int`         — accepted and ignored (belonged to the superseded
                               Taylor trig PCE; see `affine_trig_pce`)
- `angle_params`  — `Dict` mapping each load unit's ID (Int) to its angle
                    range `(v_min_deg::Float64, v_max_deg::Float64)` in degrees
                    at harmonic `h`.

# Returns
`NamedTuple` with fields:
- `I_star_independent`, `I_star_shared`                — optimal I_h [pu]
- `C_scalar_independent`, `C_scalar_shared`            — Cantelli coefficient vectors
- `binding_bus_independent`, `binding_bus_shared`      — binding bus IDs

!!! note "Voltage reference"
    The IHD limits come from `bus_voltage_limits(data, h)`, which uses
    `bus["vm"]` from `data` — the **nominal** fundamental voltage on raw parsed
    data. The JuMP `sHHC_SOC` solver uses the **HPF-solved** fundamental, so its
    `I_h` is typically a few percent lower on networks with voltage drop. To
    align this analysis with the solver, pre-populate `data` bus `vm` with HPF
    results before calling. See [`bus_voltage_limits`](@ref).
"""
function compare_shared_vs_independent(data        ::Dict,
                                        h           ::Int,
                                        epsilon     ::Float64,
                                        beta_params ::NamedTuple,
                                        angle_params::Dict;
                                        bound       ::Symbol = :cantelli)
    α̂       = beta_params.alpha_shape
    β̂       = beta_params.beta_shape
    δ       = beta_params.delta
    n_terms = hasproperty(beta_params, :n_terms) ? beta_params.n_terms : 4

    load_entries = sort([(id  = parse(Int, k),
                          bus = v["load_bus"])
                         for (k, v) in data["load"]], by = x -> x.id)
    unit_ids     = [l.id  for l in load_entries]
    load_bus_ids = [l.bus for l in load_entries]
    n_loads      = length(unit_ids)

    Z_h = harmonic_transfer_impedance(data, h, load_bus_ids)

    active_bus_ids = _active_bus_ids_lp(data, h)

    volt_lim_sq = bus_voltage_limits(data, h)

    # (a) Independent germs — N = n_loads germ dimensions
    pb_ind = PCEBasis(n_loads, α̂, β̂, δ)

    # Each unit occupies its own germ dimension (position in units vector = germ index)
    units_ind = [(id      = unit_ids[u],
                  h_angles= Dict(h => angle_params[unit_ids[u]]))
                 for u in 1:n_loads]

    cos_c_ind, sin_c_ind = angle_pce_coefficients(pb_ind, units_ind, [h];
                                                  n_terms = n_terms)

    C_re_ind, C_im_ind = compute_C_coefficients(
        Z_h, cos_c_ind, sin_c_ind, unit_ids, pb_ind, h
    )
    alpha_ind    = compute_alpha(C_re_ind, C_im_ind, pb_ind)
    Cs_ind       = compute_C_scalar(alpha_ind, pb_ind, epsilon; bound = bound)
    I_ind, b_ind = lp_solution(Cs_ind, volt_lim_sq, active_bus_ids)

    # (b) Shared germ — N = 1 germ dimension, representative angle distribution
    pb_sh = PCEBasis(1, α̂, β̂, δ)

    # Use the first unit's angle bounds as the representative shared distribution
    first_uid = unit_ids[1]
    ref_unit  = [(id      = first_uid,
                  h_angles= Dict(h => angle_params[first_uid]))]

    cos_c_ref, sin_c_ref = angle_pce_coefficients(pb_sh, ref_unit, [h];
                                                  n_terms = n_terms)
    shared_cos = cos_c_ref[(first_uid, h)]
    shared_sin = sin_c_ref[(first_uid, h)]

    # Broadcast shared PCE vectors to all unit IDs
    cos_c_sh = Dict((uid, h) => shared_cos for uid in unit_ids)
    sin_c_sh = Dict((uid, h) => shared_sin for uid in unit_ids)

    C_re_sh, C_im_sh = compute_C_coefficients(
        Z_h, cos_c_sh, sin_c_sh, unit_ids, pb_sh, h
    )
    alpha_sh   = compute_alpha(C_re_sh, C_im_sh, pb_sh)
    Cs_sh      = compute_C_scalar(alpha_sh, pb_sh, epsilon; bound = bound)
    I_sh, b_sh = lp_solution(Cs_sh, volt_lim_sq, active_bus_ids)

    return (
        I_star_independent       = I_ind,
        I_star_shared            = I_sh,
        C_scalar_independent     = Cs_ind,
        C_scalar_shared          = Cs_sh,
        binding_bus_independent  = b_ind,
        binding_bus_shared       = b_sh,
    )
end
