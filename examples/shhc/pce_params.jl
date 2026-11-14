################################################################################
# HarmonicPowerModels.jl                                                       #
# Shared PCE parameters for the SHHC paper examples                            #
# Beta distribution fitted to 6-pulse LCC converter harmonic phase angles      #
# See: §II of the paper, Table III                                              #
################################################################################

const PCE_PARAMS = Dict(
    # Beta distribution shape parameters (shared across all harmonics)
    "alpha_hat" => 4.0414,
    "beta_hat"  => 4.7762,

    # PCE degree (fixed at δ=2 throughout the paper)
    "delta"     => 2,

    # Taylor series terms for trig approximation (n=4 is adequate, see Table II)
    "n_terms"   => 4,

    # Primary risk level: ε=0.05 consistent with EN50160 (5% violation allowed)
    "epsilon"   => 0.05,

    # Harmonic-specific location (v_min) and scale (v_max) in radians.
    # Derived from §II Monte Carlo analysis of a representative 6-pulse LCC unit.
    # Keys are harmonic orders (Int).
    "harmonics" => Dict(
        5  => (v_min = deg2rad(  4.75), v_max = deg2rad( 40.05)),  # range = 35.30°
        7  => (v_min = deg2rad(-22.77), v_max = deg2rad( 27.51)),  # range = 50.28°
        11 => (v_min = deg2rad( 79.69), v_max = deg2rad(164.01)),  # range = 84.32°
        13 => (v_min = deg2rad( 49.38), v_max = deg2rad(147.67)),  # range = 98.29°
    ),
)

# Phase-controlled (thyristor) bridge, §II-D, Table IV: firing delay, L_ac and L_dc
# sampled jointly, i.e., an independent firing delay per unit. Selected with
# ENV["SHHC_UNIT"] = "thyristor".
const PCE_PARAMS_THYRISTOR = Dict(
    "alpha_hat" => 4.139,
    "beta_hat"  => 3.659,
    "harmonics" => Dict(
        5  => (v_min = deg2rad(-138.23), v_max = deg2rad(-10.87)),  # range = 127.37°
        7  => (v_min = deg2rad(-265.07), v_max = deg2rad(-48.90)),  # range = 216.17°
        11 => (v_min = deg2rad(-240.99), v_max = deg2rad( 44.65)),  # range = 285.64°
        13 => (v_min = deg2rad(-337.52), v_max = deg2rad( 20.25)),  # range = 357.76°
    ),
)
# Extension to the 25th harmonic (Jinwen Liu, supplementary tables of the paper,
# 2026-11-11): the shapes are refitted jointly over h = 5, ..., 25, so the 5th to
# 13th harmonic ranges change slightly as well. Selected with ENV["SHHC_HMAX"] = "25".
const SHHC_H25 = get(ENV, "SHHC_HMAX", "13") == "25"
if SHHC_H25
    merge!(PCE_PARAMS, Dict(
        "alpha_hat" => 3.731,
        "beta_hat"  => 4.439,
        "harmonics" => Dict(
            5  => (v_min = deg2rad(   4.80), v_max = deg2rad( 40.90)),
            7  => (v_min = deg2rad( -22.52), v_max = deg2rad( 27.18)),
            11 => (v_min = deg2rad(  77.87), v_max = deg2rad(163.96)),
            13 => (v_min = deg2rad(  48.98), v_max = deg2rad(148.35)),
            17 => (v_min = deg2rad(-239.51), v_max = deg2rad(-76.70)),
            19 => (v_min = deg2rad(  80.08), v_max = deg2rad(266.23)),
            23 => (v_min = deg2rad(-209.79), v_max = deg2rad( 34.31)),
            25 => (v_min = deg2rad(-246.50), v_max = deg2rad( 15.56)),
        )))
    merge!(PCE_PARAMS_THYRISTOR, Dict(
        "alpha_hat" => 4.1612,
        "beta_hat"  => 3.5541,
        "harmonics" => Dict(
            5  => (v_min = deg2rad(-143.07), v_max = deg2rad( -9.30)),
            7  => (v_min = deg2rad(-263.84), v_max = deg2rad(-52.49)),
            11 => (v_min = deg2rad(-248.03), v_max = deg2rad( 46.71)),
            13 => (v_min = deg2rad(-353.10), v_max = deg2rad(  5.22)),
            17 => (v_min = deg2rad(-346.15), v_max = deg2rad(104.27)),
            19 => (v_min = deg2rad(-447.14), v_max = deg2rad( 62.43)),
            23 => (v_min = deg2rad(-442.00), v_max = deg2rad(161.32)),
            25 => (v_min = deg2rad(-180.49), v_max = deg2rad(479.68)),
        )))
end
get(ENV, "SHHC_UNIT", "diode") == "thyristor" && merge!(PCE_PARAMS, PCE_PARAMS_THYRISTOR)

# Risk levels for the sensitivity study in §IV.B
const EPSILON_SENSITIVITY = [0.05, 0.10, 0.20]
