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
        form = NLP_DHHC

        # define the set of considered harmonics
        H=[1, 3, 5, 7, 9, 13]

        # build harmonic data
        hdata = HPM.replicate(data, H=H)

        # solve HC problem
        results_hhc = HPM.solve_hhc(hdata, form, solver)

        @testset "General Tests" begin
            # Solved to optimality
            @test results_hhc["termination_status"] == LOCALLY_SOLVED
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

        @testset "Transformer Phase Shifts" begin
            # The phase shift for the harmonic voltage over a transformer, i.e.,
            # between the excitation Eₓₕ and internal voltage of the second 
            # winding V₂ₕ in a balanced (optimal) power flow should relate as 
            # follows:
            # - 'positive sequence' harmonics, i.e., h = 1 + 3n = 1, 4, 7...
            #       θⱽ - θᵉ ≈ vg * π/6
            # - 'negative sequence' harmonics, i.e., h = 2 + 3n = 2, 5, 8...
            #       θⱽ - θᵉ ≈ -vg * π/6
            # - 'zero sequence' harmonics, i.e., h = 3 + 3n = 3, 6, 9...
            #       θⱽ ≈ θᵉ
            for nh ∈ H, (nx,xfmr) in data["xfmr"]
                vg  = parse(Int, xfmr["vg"][3:end])
                
                sol_xfmr = results_hhc["solution"]["nw"]["$nh"]["xfmr"][nx]

                θe = angle(sol_xfmr["ert"] + im * sol_xfmr["eit"])
                θv = angle(sol_xfmr["vrt_to"] + im * sol_xfmr["vit_to"])

                if HPM.is_pos_sequence(nh)
                    @test isapprox(rem(θv - θe - vg * π/6, 2π), 0.0, atol=1e-6)
                elseif HPM.is_neg_sequence(nh)
                    @test isapprox(rem(θv - θe + vg * π/6, 2π), 0.0, atol=1e-6)
                elseif HPM.is_zero_sequence(nh)
                    @test isapprox(rem(θv - θe, 2π), 0.0, atol=1e-6)
                end
            end
        end
    end
end