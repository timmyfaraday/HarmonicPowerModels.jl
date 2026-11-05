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
get(ENV, "SHHC_UNIT", "diode") == "thyristor" && merge!(PCE_PARAMS, PCE_PARAMS_THYRISTOR)

# Risk levels for the sensitivity study in §IV.B
const EPSILON_SENSITIVITY = [0.05, 0.10, 0.20]
