################################################################################
# HarmonicPowerModels.jl — converter types and type mixes for the IEEE30         #
# recombination study (feature/shhc-recombination, Phase 2)                     #
#                                                                               #
# Every harmonic angle here is relative to the unit's own fundamental bus       #
# voltage (local angle reference). All types keep the shared Beta(α̂, β̂) shape   #
# and one germ per unit, so the number of random variables is unchanged.        #
#                                                                               #
# Types (assumptions — the paper's §II data covers type A only):                #
#   A      the §II 6-pulse line-commutated converter: Table III ranges.         #
#   B(αf)  the same converter with phase control at firing delay αf. Delaying   #
#          the firing delays the whole current waveform by αf, which rotates    #
#          the h-th harmonic phasor of phase a by −h·αf; the dispersion (range  #
#          width) is kept. Exact for the idealised rectangular current, an      #
#          approximation once commutation overlap changes with αf.              #
#   B(αc, w) type B whose firing delay is not fixed but ranges over          #
#          [αc − w/2, αc + w/2] with the operating point. Approximation that    #
#          keeps one germ per unit: the mean shifts by −h·αc and the support    #
#          widens by h·w (the delay range seen at harmonic h), with the shared  #
#          Beta shape over the widened support.                                 #
#   C      a PWM active front end. Its low-order harmonics are not set by       #
#          commutation and their angle follows the control and background       #
#          distortion, so it is modelled as spread over the full circle:        #
#          [E_A − 180°, E_A + 180°] around the type-A mean.                     #
#                                                                               #
# Mixes assign types to the harmonic units in load-index order, cycling, so      #
# every region of the network holds every type (units at comparable electrical  #
# positions differ in technology, not in location):                             #
#   M0  all A                                   (Phase 1 reference)             #
#   M1  A, B(30°) alternating                                                   #
#   M2  A, B(15°), B(30°), B(45°) cycling        (spread of firing delays)      #
#   M3  A, C alternating                                                        #
# Author: Tom Van Acker                                                          #
################################################################################

const TYPE_A = Dict(h => (PCE_PARAMS["harmonics"][h].v_min, PCE_PARAMS["harmonics"][h].v_max)
                    for h in (5, 7, 11, 13))

"Angle range [rad] of converter `type` at harmonic `h`: :A, (:B, αf_deg), (:B, αc_deg, w_deg) or :C."
function type_range(type, h)
    lo, hi = TYPE_A[h]
    type == :A && return (lo, hi)
    if type isa Tuple && type[1] == :B && length(type) == 2
        δ = -h * deg2rad(type[2])
        return (lo + δ, hi + δ)
    end
    if type isa Tuple && type[1] == :B && length(type) == 3
        Eω   = PCE_PARAMS["alpha_hat"] / (PCE_PARAMS["alpha_hat"] + PCE_PARAMS["beta_hat"])
        mean = lo + (hi - lo) * Eω - h * deg2rad(type[2])
        W    = (hi - lo) + h * deg2rad(type[3])
        return (mean - W * Eω, mean + W * (1 - Eω))
    end
    if type == :C
        Eω = PCE_PARAMS["alpha_hat"] / (PCE_PARAMS["alpha_hat"] + PCE_PARAMS["beta_hat"])
        mean = lo + (hi - lo) * Eω
        # support of 2π placed so that its Beta mean equals the type-A mean
        return (mean - 2π * Eω, mean + 2π * (1 - Eω))
    end
    error("unknown converter type $type")
end

const MIXES = Dict(
    "M0" => [:A],
    "M1" => [:A, (:B, 30.0)],
    "M2" => [:A, (:B, 15.0), (:B, 30.0), (:B, 45.0)],
    "M3" => [:A, :C],
)

"Per-harmonic, per-unit angle ranges for `mix` over the sorted `units`."
function mix_ranges(mix::String, units)
    cyc = MIXES[mix]
    Dict{Int,Any}(h => Dict{Int,Tuple{Float64,Float64}}(
                      d => type_range(cyc[mod1(k, length(cyc))], h) for (k, d) in enumerate(units))
                  for h in (5, 7, 11, 13))
end

"""
Per-unit angle ranges with a share `s` ∈ [0, 1] of the units of type `typeB`
and the rest of type A. Type-B units are spread evenly over the load order:
unit k (1-based) is type B when ⌊k·s⌋ > ⌊(k−1)·s⌋.
"""
function share_ranges(s::Real, typeB, units)
    isB(k) = floor(Int, k * s + 1e-9) > floor(Int, (k - 1) * s + 1e-9)
    Dict{Int,Any}(h => Dict{Int,Tuple{Float64,Float64}}(
                      d => type_range(isB(k) ? typeB : :A, h) for (k, d) in enumerate(units))
                  for h in (5, 7, 11, 13))
end
