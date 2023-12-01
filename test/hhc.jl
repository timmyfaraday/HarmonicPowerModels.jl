################################################################################
#  Copyright 2023, Tom Van Acker                                               #
################################################################################
# HarmonicPowerModels.jl                                                       #
# An extention package of PowerModels(Distribution).jl for Harmonics           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################

@testset "Harmonic Hosting Capacity" begin
    @testset "Industrial Network" begin

        # read-in data 
        path = joinpath(HPM.BASE_DIR,"test/data/matpower/industrial_network_hhc.m")
        data = PMs.parse_file(path)

        # set the formulation
        form = PMs.IVRPowerModel

        # define the set of considered harmonics
        H=[1, 3, 5, 7, 9, 13]

        # build harmonic data
        hdata = HPM.replicate(data, H=H)
        
        # Define ref_angle and angle range 
        for nh in H
            for (l, load) in hdata["nw"]["$nh"]["load"]
                load["reference_harmonic_angle"] = 0 # pi / 4 # rad
                load["harmonic_angle_range"] = 0 # pi / 10 # rad, symmetric around reference
            end
        end

        # solve HC problem
        results_hhc = HPM.solve_hhc(hdata, form, solver)
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

        # @testset "Transformer Phase Shifts" begin
            
        #     for nh ∈ H, (nx,xfmr) in data["xfmr"]
        #         vg  = parse(Int, xfmr["vg"][3:end])
                


        #         sol_xfmr = results_hhc["solution"]["nw"]["$nh"]["xfmr"][nx]

        #         θe = rad2deg(angle(sol_xfmr["ert"] + im * sol_xfmr["eit"]))
        #         θv = rad2deg(angle(sol_xfmr["vrt_to"] + im * sol_xfmr["vit_to"]))

        #         @test θv - θa ⪅
        #     end
        # end
    end
end