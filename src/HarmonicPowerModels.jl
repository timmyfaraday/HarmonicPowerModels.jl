################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Authors: Tom Van Acker, Frederik Geth, Hakan Ergun                           #
################################################################################
# Changelog:                                                                   #
# v0.2.0 - reviewed TVA                                                        #
# v0.2.1 - reviewed TVA                                                        #
################################################################################

module HarmonicPowerModels

    # using pkgs 
    using ProgressMeter

    # import pkgs
    import JuMP
    import MathOptInterface

    import PowerModels
    import InfrastructureModels

    import SignalDecomposition
    import Interpolations
    import Distributions
    import Random
    import Statistics

    # import function to overwrite
    import InfrastructureModels: replicate

    # pkg constants 
    const _HPM = HarmonicPowerModels

    const _MOI = MathOptInterface

    const _PMs = PowerModels
    const _IMs = InfrastructureModels
    
    const _SDC = SignalDecomposition
    const _INT = Interpolations

    # const 
    const freq = 50.0

    # funct
    fundamental(pm) = 1
    sorted_nw_ids(pm) = sort(collect(_PMs.nw_ids(pm)))

    # paths
    const BASE_DIR = dirname(@__DIR__)

    # include
    include("core/base.jl")
    include("core/pce_basis.jl")    # unified PolyChaos-based PCE workspace
    include("core/types.jl")
    include("core/constraint_template.jl")
    include("core/data.jl")
    include("core/variable.jl")

    include("form/iv.jl")
    include("form/iv_pce.jl")

    include("prob/hopf.jl")
    include("prob/hpf.jl")
    include("prob/hhc.jl")
    include("prob/shhc.jl")

    include("util/imp.jl")
    include("util/init.jl")
    include("util/ref.jl")
    include("util/sol.jl")
    include("util/xfmr_magn.jl")
    include("util/harmonic_admittance.jl")
    include("core/lp_reduction.jl")

    # export
    export BASE_DIR
    export dHHC_NLP, dHHC_SOC
    export sHHC_QCQP, sHHC_SOC

    export replicate
    export solve_hpf, solve_hopf, solve_hhc, solve_shhc
    export update_hdata_with_pce_data!, update_hdata_with_angle_offsets!

    export calculate_pos_seq_harmonic_impedance
    export build_harmonic_Ybus, harmonic_transfer_impedance, bus_voltage_limits
    export reduced_harmonic_impedance, constructive_worst_case_angles, recentre_angle_ranges
    export empirical_lambdas, solve_shhc_calibrated
    export compute_C_coefficients, compute_alpha, compute_C_scalar
    export lp_solution, compare_shared_vs_independent
    
end
