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
# v0.3.0 - adapted for extend graph representation                             #
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
    import SparseArrays

    import Memento
    
    # pkg constants 
    const _HPM = HarmonicPowerModels

    const _MOI = MathOptInterface

    const _PMs = PowerModels
    const _IMs = InfrastructureModels
    
    const _SDC = SignalDecomposition
    const _INT = Interpolations
    const _SPA = SparseArrays

    const _MEM = Memento

    # const 
    const freq = 50.0

    # funct
    fundamental(pm) = 1
    sorted_nw_ids(pm) = sort(collect(_PMs.nw_ids(pm)))

    # paths
    const BASE_DIR = dirname(@__DIR__)

    # include
    ## core
    include("core/base.jl")
    include("core/types.jl")
    include("core/data.jl")
    include("core/objective.jl")
    include("core/solution.jl")

    ## cmp - bus
    include("cmp/bus/bus.jl")
    include("cmp/bus/bus_ref.jl")
    include("cmp/bus/bus_clean.jl")
    ## cmp - edge 
    include("cmp/edge/branch.jl")
    include("cmp/edge/xfmr.jl")
    ## cmp - unit
    include("cmp/unit/filter.jl")
    include("cmp/unit/gen.jl")
    include("cmp/unit/hload.jl")
    include("cmp/unit/hsrc.jl")
    include("cmp/unit/shunt.jl")

    ## prob
    include("prob/hopf.jl")
    include("prob/hpf.jl")
    include("prob/hhc.jl")

    ## util 
    include("util/ihd.jl")
    include("util/imp.jl") 
    include("util/ref.jl")
    include("util/thd.jl")
    include("util/xfmr_magn.jl")

    # export
    export BASE_DIR
    export HarmonicPowerModel, dHHCPowerModel

    export build_hdata_from_matpower_file
    export solve_hpf, solve_hopf, solve_hhc

    export calculate_pos_seq_harmonic_impedance
end
