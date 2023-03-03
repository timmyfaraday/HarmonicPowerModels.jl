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
    import PowerModelsDistribution
    import SignalDecomposition
    import InfrastructureModels
    import Interpolations
    import Memento

    # import types
    import PowerModels: AbstractPowerModel, AbstractIVRModel
    import PowerModels: ids, ref, var, con, sol, nw_ids, nws
    import InfrastructureModels: replicate, sol_component_value_edge

    # pkg constants 
    const _PMs = PowerModels
    const _PMD = PowerModelsDistribution
    const _HPM = HarmonicPowerModels
    const _SDC = SignalDecomposition
    const _IMs = InfrastructureModels
    const _INT = Interpolations

    function __init__()
        global _LOGGER = Memento.getlogger(PowerModels)
    end

    # const 
    const freq = 50.0
    const nw_id_default = 1

    # funct
    sorted_nw_ids(pm) = sort(collect(_PMs.nw_ids(pm)))

    # paths
    const BASE_DIR = dirname(@__DIR__)

    # include
    include("core/constraint_template.jl")
    include("core/data.jl")
    include("core/variable.jl")

    include("form/ivr.jl")

    include("prob/hopf.jl")
    include("prob/hpf.jl")

    include("util/sol.jl")
    include("util/xfmr_magn.jl")
    include("util/xfmr.jl")

    # export
    export BASE_DIR

    export replicate

    export solve_hopf, solve_hpf

end
