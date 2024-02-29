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
################################################################################

module HarmonicPowerModels

    # using pkgs 
    using ProgressMeter

    # import pkgs
    import JuMP
    import PowerModels
    import InfrastructureModels

    import SignalDecomposition
    import Interpolations
    
    # import function to overwrite
    import InfrastructureModels: replicate

    # pkg constants 
    const _HPM = HarmonicPowerModels

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
    include("core/types.jl")
    include("core/constraint_template.jl")
    include("core/data.jl")
    include("core/variable.jl")

    include("form/iv.jl")

    include("prob/hopf.jl")
    include("prob/hpf.jl")
    include("prob/hhc.jl")

    include("util/init.jl")
    include("util/ref.jl")
    include("util/sol.jl")
    include("util/xfmr_magn.jl")

    # export
    export BASE_DIR
    export dHHC_NLP, dHHC_SOC

    export replicate
    export solve_hpf, solve_hopf, solve_hhc 

end
