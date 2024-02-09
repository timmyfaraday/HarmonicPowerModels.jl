################################################################################
#  Copyright 2023, Frederik Geth, Tom Van Acker                                #
################################################################################
# HarmonicPowerModels.jl                                                       #
# An extention package of PowerModels(Distribution).jl for Harmonics           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
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
    
    # import types
    import PowerModels: AbstractPowerModel, AbstractIVRModel
    import PowerModels: ids, ref, var, con, sol
    import InfrastructureModels: replicate, sol_component_value_edge

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

    include("util/sol.jl")
    include("util/xfmr_magn.jl")
    include("util/xfmr.jl")
    include("util/pf.jl")

    # export
    export BASE_DIR

    export dHHC_NLP, dHHC_SOC

    export replicate

    export solve_hpf, solve_hopf, solve_hhc 

end
