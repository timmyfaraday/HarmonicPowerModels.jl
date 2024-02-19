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

    @testset "NLP - Industrial Network" begin

        # read-in data 
        path = joinpath(HPM.BASE_DIR,"test/data/matpower/industrial_network_hhc.m")
        data = PMs.parse_file(path)

        # set the formulation
        form = dHHC_NLP

        # define the set of considered harmonics
        H   = [1, 3, 5, 7, 9, 13]
        H⁰  = H[(H.%3).==0]

        # build harmonic data
        hdata = HPM.replicate(data, H=H)

        # solve HHC problem
        results_hhc = HPM.solve_hhc(hdata, form, solver_nlp)

        @testset "Feasibility and Objective" begin
            # solved to optimality
            @test results_hhc["termination_status"] == LOCALLY_SOLVED
            # objective value depending on fairness principle
            if !haskey(hdata, "principle")
                @test isapprox(results_hhc["objective"], 0.160236; atol = 1e-4)
            elseif hdata["principle"] == "equality"
                @test isapprox(results_hhc["objective"], 0.019055; atol = 1e-4)
            end
        end

        @testset "Root Mean Square Voltage" begin 
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
        
        @testset "Total Harmonic Voltage Distortion" begin
            # THDᵢ = √(∑ₕ(|Uᵢₕ|²) / |Uᵢ₁|²) ≤ THDmaxᵢ, ∀ i ∈ I
            for (nb,bus) ∈ data["bus"]
                thdmax  = bus["thdmax"]

                vm_fund = results_hhc["solution"]["nw"]["1"]["bus"][nb]["vm"]
                vm_harm = [results_hhc["solution"]["nw"]["$nh"]["bus"][nb]["vm"] 
                                for nh ∈ H if nh ≠ 1]

                @test sqrt(sum(vm_harm.^2) / vm_fund^2) ⪅ thdmax
            end
        end
        
        @testset "Individual Harmonic Voltage Distortion" begin
            # IHDᵢₕ = √(|Uᵢₕ|² / |Uᵢ₁|²) ≤ IHDmaxᵢₕ, ∀ i ∈ I, h ∈ H
            for nh ∈ H, (nb,bus) in hdata["nw"]["$nh"]["bus"] if nh ≠ 1
                ihdmax  = bus["ihdmax"]

                vm_fund = results_hhc["solution"]["nw"]["1"]["bus"][nb]["vm"]
                vm_harm = results_hhc["solution"]["nw"]["$nh"]["bus"][nb]["vm"]

                @test sqrt(vm_harm^2 / vm_fund^2) ⪅ ihdmax
            end end
        end
        
        @testset "Zero-Sequence Current Blocking" begin 
            # |I|ₓᵢⱼₕ = 0, ∀ xij ∈ Tˣ⁻ⁿᵉ ∪ Tˣ⁻ᵈᵉˡᵗᵃ, h ∈ H⁰
            # |Iˢ|ₓᵢⱼₕ = 0, ∀ xij ∈ Tˣ⁻ⁿᵉ, h ∈ H⁰
            for nh ∈ H⁰, (nx, xfmr) in hdata["nw"]["$nh"]["xfmr"]
                cnf1    = xfmr["cnf1"]
                cnf2    = xfmr["cnf2"]

                gnd1    = xfmr["gnd1"]
                gnd2    = xfmr["gnd2"]

                ct_fr   = results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["crt_fr"] +
                          results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["cit_fr"] * im
                ct_to   = results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["crt_to"] + 
                          results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["cit_to"] * im 

                cst_fr  = results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["csrt_fr"] +
                          results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["csit_fr"] * im
                cst_to  = results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["csrt_to"] + 
                          results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["csit_to"] * im

                if cnf1 == 'D' 
                    @test abs(ct_fr) ≈ 0.0
                end
                if cnf1 ≠ 'D' && gnd1 == 0 
                    @test abs(ct_fr) ≈ 0.0
                    @test abs(cst_fr) ≈ 0.0
                end
                if cnf2 == 'D' 
                    @test abs(ct_to) ≈ 0.0
                end
                if cnf2 ≠ 'D' && gnd2 == 0 
                    @test abs(ct_to) ≈ 0.0
                    @test abs(cst_to) ≈ 0.0
                end
            end
        end

        @testset "Branch Currents" begin
            # Īᵦᵢⱼₕ = Ūᵢ * yˢʰᵦᵢⱼₕ ± Īˢᵦᵢⱼₕ, ∀ βij in Tᵇ, h in H
            # Īˢᵦᵢⱼₕ = (Ūᵢₕ - Ūⱼₕ) * zᵦ, ∀ βij in Tᵇ⁻ᶠʳ, h in H
            for nh ∈ H, (nb, branch) in hdata["nw"]["$nh"]["branch"]
                f_bus   = branch["f_bus"]
                t_bus   = branch["t_bus"]

                y_fr    = branch["g_fr"] + im * branch["b_fr"]
                y_to    = branch["g_to"] + im * branch["b_to"]

                z       = branch["br_r"] + im * branch["br_x"]

                v_fr    = results_hhc["solution"]["nw"]["$nh"]["bus"]["$f_bus"]["vr"] +
                          results_hhc["solution"]["nw"]["$nh"]["bus"]["$f_bus"]["vi"] * im
                v_to    = results_hhc["solution"]["nw"]["$nh"]["bus"]["$t_bus"]["vr"] + 
                          results_hhc["solution"]["nw"]["$nh"]["bus"]["$t_bus"]["vi"] * im
                
                cs      = results_hhc["solution"]["nw"]["$nh"]["branch"]["$nb"]["csr_fr"] +
                          results_hhc["solution"]["nw"]["$nh"]["branch"]["$nb"]["csi_fr"] * im
                c_fr    = results_hhc["solution"]["nw"]["$nh"]["branch"]["$nb"]["cr_fr"] + 
                          results_hhc["solution"]["nw"]["$nh"]["branch"]["$nb"]["ci_fr"] * im
                c_to    = results_hhc["solution"]["nw"]["$nh"]["branch"]["$nb"]["cr_to"] +
                          results_hhc["solution"]["nw"]["$nh"]["branch"]["$nb"]["ci_to"] * im

                # branch total current 
                @test c_fr ≈ v_fr * y_fr + cs
                @test c_to ≈ v_to * y_to - cs
                # Ohm's law 
                @test v_fr - v_to ≈ z * cs
            end
        end

        @testset "Transformer Voltages and Currents" begin
            # Ēₓₕ = tᵛᵍₓₕ * Vₓⱼᵢₕ, ∀ xji ∈ Tˣ⁻ᵗᵒ, h ∈ H 
            # Īᵐₓₕ + Ēₓₕ = Īˢₓᵢⱼₕ + tᵛᵍₓₕ * Īₓⱼᵢₕ, ∀ xij ∈ Tˣ⁻ᶠʳ, h ∈ H 
            for nh ∈ H, (nx, xfmr) in hdata["nw"]["$nh"]["xfmr"]
                f_bus   = xfmr["f_bus"]
                t_bus   = xfmr["t_bus"]

                t_vg    = xfmr["tr"] + im * xfmr["ti"]

                v_fr    = results_hhc["solution"]["nw"]["$nh"]["bus"]["$f_bus"]["vr"] +
                          results_hhc["solution"]["nw"]["$nh"]["bus"]["$f_bus"]["vi"] * im
                v_to    = results_hhc["solution"]["nw"]["$nh"]["bus"]["$t_bus"]["vr"] + 
                          results_hhc["solution"]["nw"]["$nh"]["bus"]["$t_bus"]["vi"] * im

                rsh     = xfmr["rsh"]
                z       = xfmr["xsc"] * im

                e       = results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["ert"] +
                          results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["eit"] * im

                vt_fr   = results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["vrt_fr"] +
                          results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["vit_fr"] * im
                vt_to   = results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["vrt_to"] +
                          results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["vit_to"] * im

                cmt     = results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["cmrt"] +
                          results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["cmit"] * im

                ct_fr   = results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["crt_fr"] +
                          results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["cit_fr"] * im
                ct_to   = results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["crt_to"] + 
                          results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["cit_to"] * im 

                cst_fr  = results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["csrt_fr"] +
                          results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["csit_fr"] * im
                cst_to  = results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["csrt_to"] + 
                          results_hhc["solution"]["nw"]["$nh"]["xfmr"]["$nx"]["csit_to"] * im

                println("h=$nh, for xfmr $nx")

                # Kirchhoff's current law
                @test cmt + e / rsh ≈ cst_fr + conj(t_vg) * cst_to # why conj?
                # Ohm's law
                @test vt_fr ≈ e + z * cst_fr
                # transformer phase shift
                @test e ≈ conj(t_vg) * vt_to # why conj?
            end
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
            if !haskey(hdata, "principle")
                @test isapprox(results_hhc_soc["objective"], 0.160236; atol = 1e-4)
            elseif hdata["principle"] == "equality"
                @test isapprox(results_hhc_soc["objective"], 0.019055; atol = 1e-4)
            end
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