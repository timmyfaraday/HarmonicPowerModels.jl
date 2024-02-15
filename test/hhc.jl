################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Authors: Tom Van Acker, Hakan Ergun                                          #
################################################################################
# Changelog:                                                                   #
################################################################################

@testset "Harmonic Hosting Capacity" begin

    @testset "Industrial Network" begin

        # read-in data 
        path = joinpath(HPM.BASE_DIR,"test/data/matpower/industrial_network_hhc.m")
        data = PMs.parse_file(path)

        # set the formulation
        form = dHHC_NLP

        # define the set of considered harmonics
        H=[1, 3, 5, 7, 9, 13]

        # build harmonic data
        hdata = HPM.replicate(data, H=H)

        # solve HHC problem
        results_hhc = HPM.solve_hhc(hdata, form, solver_nlp)

        @testset "Feasibility" begin
            # Solved to optimality
            @test results_hhc["termination_status"] == LOCALLY_SOLVED
            @test isapprox(results_hhc["objective"], 0.160237; atol = 1e-4)
        end

        @testset "Root Mean Square" begin 
            # Uminᵢ ≤ RMSᵢ = √(∑ₕ(|Uᵢₕ|²)) ≤ Umaxᵢ, ∀ i ∈ I
            for (nb, bus) ∈ data["bus"]
                vmin    = bus["vmin"]
                vmax    = bus["vmax"]
                
                vm      = [results_hhc["solution"]["nw"]["$nh"]["bus"][nb]["vm"] 
                                for nh ∈ H]

                @test vmin ⪅ sqrt(sum(vm.^2))
                @test sqrt(sum(vm.^2)) ⪅ vmax
            end
        end
        
        @testset "Total Harmonic Distortion" begin
            # THDᵢ = √(∑ₕ(|Uᵢₕ|²) / |Uᵢ₁|²) ≤ THDmaxᵢ, ∀ i ∈ I
            for (nb,bus) ∈ data["bus"]
                thdmax  = bus["thdmax"]

                vm_fund = results_hhc["solution"]["nw"]["1"]["bus"][nb]["vm"]
                vm_harm = [results_hhc["solution"]["nw"]["$nh"]["bus"][nb]["vm"] 
                                for nh ∈ H if nh ≠ 1]

                @test sqrt(sum(vm_harm.^2) / vm_fund^2) ⪅ thdmax
            end
        end
        
        @testset "Individual Harmonic Distortion" begin
            # IHDᵢₕ = √(|Uᵢₕ|² / |Uᵢ₁|²) ≤ IHDmaxᵢₕ, ∀ i ∈ I, h ∈ H
            for nh ∈ H, (nb,bus) in hdata["nw"]["$nh"]["bus"] if nh ≠ 1
                ihdmax  = bus["ihdmax"]

                vm_fund = results_hhc["solution"]["nw"]["1"]["bus"][nb]["vm"]
                vm_harm = results_hhc["solution"]["nw"]["$nh"]["bus"][nb]["vm"]

                @test sqrt(vm_harm^2 / vm_fund^2) ⪅ ihdmax
            end end
        end    
    end

    @testset "Industrial Network NLP vs SOC" begin
        # read-in data 
        path = joinpath(HPM.BASE_DIR,"test/data/matpower/industrial_network_hhc.m")
        data = PMs.parse_file(path)

        # set the formulation
        form = dHHC_SOC

        # define the set of considered harmonics
        H=[1, 3, 5, 7, 9, 13]

        # build harmonic data
        hdata = HPM.replicate(data, H=H)
        for (n, nw) in hdata["nw"]
            for (b, bus) in nw["bus"]
                bus["angle_range"] = 0.0
            end
        end

        # solve HHC problem
        results_hhc_soc = HPM.solve_hhc(hdata, form, solver_soc, solver_nlp)

        @testset "Feasibility" begin
            @test results_hhc_soc["termination_status"] == ALMOST_OPTIMAL
            @test isapprox(results_hhc_soc["objective"], 0.160236; atol = 1e-4)
        end

        @testset "Root Mean Square" begin 
            # Uminᵢ ≤ RMSᵢ = √(∑ₕ(|Uᵢₕ|²)) ≤ Umaxᵢ, ∀ i ∈ I
            for (nb, bus) ∈ data["bus"]
                vmin    = bus["vmin"]
                vmax    = bus["vmax"]
                
                vm      = vcat(hdata["nw"]["1"]["bus"][nb]["vm"], [results_hhc_soc["solution"]["nw"]["$nh"]["bus"][nb]["vm"] for nh ∈ setdiff(H,1)])

                @test vmin ⪅ sqrt(sum(vm.^2))
                @test sqrt(sum(vm.^2)) ⪅ vmax
            end
        end

        @testset "Total Harmonic Distortion" begin
            # THDᵢ = √(∑ₕ(|Uᵢₕ|²) / |Uᵢ₁|²) ≤ THDmaxᵢ, ∀ i ∈ I
            for (nb,bus) ∈ data["bus"]
                thdmax  = bus["thdmax"]

                vm_fund = hdata["nw"]["1"]["bus"][nb]["vm"]
                vm_harm = [results_hhc_soc["solution"]["nw"]["$nh"]["bus"][nb]["vm"] 
                                for nh ∈ H if nh ≠ 1]

                @test sqrt(sum(vm_harm.^2) / vm_fund^2) ⪅ thdmax
            end
        end

        @testset "Individual Harmonic Distortion" begin
            # IHDᵢₕ = √(|Uᵢₕ|² / |Uᵢ₁|²) ≤ IHDmaxᵢₕ, ∀ i ∈ I, h ∈ H / 1
            for nh ∈ setdiff(H,1), (nb,bus) in hdata["nw"]["$nh"]["bus"]
                ihdmax  = bus["ihdmax"]

                vm_fund = hdata["nw"]["1"]["bus"][nb]["vm"]
                vm_harm = results_hhc_soc["solution"]["nw"]["$nh"]["bus"][nb]["vm"]

                @test sqrt(vm_harm^2 / vm_fund^2) ⪅ ihdmax
            end
        end
    end
end