    ################################################################################

    # HarmonicPowerModels.jl                                                       #
    # Extension package of PowerModels.jl for Steady-State Power System            #
    # Optimization with Power Harmonics.                                           #
    # See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
    ################################################################################
    # Example considering the harmonic hosting capacity of an industrial power     #
    # system taken from: Exact Lower Bound on Equitable Harmonic Hosting Capacity  #
    # by T. Van Acker and H. Ergun, pg. 6, § III.A.                                #
    # ---------------------------------------------------------------------------- #
    # This numerical illustration considers a section of a real-world industrial   #
    # phase-balanced three-phase power system across five voltage levels, with     #
    # four transformers, three cables, and four harmonic units.                    #
    # On top of the fundamental harmonic, the considered harmonic set is           # 
    # H ∈ {3, 5, 7}, i.e., one harmonic from each sequence subset. The aim of this # 
    # numerical illustration is to                                                 #
    #   1) validate the proposed formulations; and                                 #
    #   2) show the equivalence between the non-linear and second-order cone       #
    #      models.                                                                 #
    # For the purpose of proving equivalency, both models are solved using Ipopt,  # 
    # where the second-order constraints are translated to their equivalent        #
    # non-convex quadratic form using the appropriate MATHOPTINTERFACE constraint  #
    # bridge.                                                                      #
    ################################################################################
    # Authors: Tom Van Acker, Hakan Ergun                                          #
    ################################################################################
    # Changelog:                                                                   #
    # v0.2.0 - reviewed TVA                                                        #
    # v0.3.0 - example template                                                    # 
    ################################################################################

    # INPUT ########################################################################
    # using pkgs
    using HarmonicPowerModels, PowerModels
    using Ipopt 
    using PrettyTables

    # pkg cte
    const PMs = PowerModels
    const HPM = HarmonicPowerModels

    # set the solver
    solver_nlp = Ipopt.Optimizer

    # read-in data 
    path = joinpath(HPM.BASE_DIR,"test/data/matpower/industrial_network_hhc.m")
    data = PMs.parse_file(path)

    # define the set of considered harmonics
    H = [1, 3, 5, 7]

    hdata = HPM.build_hdata_from_matpower_file(data, H=H, prob=:hhc)
    result_opf = HPM.solve_hopf(hdata, HarmonicPowerModel, solver_nlp)
    result_pf = HPM.solve_hpf(hdata, HarmonicPowerModel, solver_nlp)
    result_hhc = HPM.solve_hhc(hdata, HarmonicPowerModel, solver_nlp)


    for (b, bus) in result_pf["solution"]["nw"]["1"]["bus"]
        println("PF: bus ", b, ", vm = ", sqrt(bus["vbr"]^2 + bus["vbi"]^2))
    end

    for (b, bus) in result_opf["solution"]["nw"]["1"]["bus"]
        println("OPF: bus ", b, " , vm = ", sqrt(bus["vbr"]^2 + bus["vbi"]^2))
    end