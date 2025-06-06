################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# This section presents a case study comparing three methods: IEEE519,         #
# IEC61000, and novel hhc optimization method, to determine the harmonic       #
# hosting capacity of a phase-balanced three-phase power system, using the     #
# IEEE 30 bus test system.                                                     #
# First, the harmonic hosting capacity is determined for a purely inductive    #
# network based on the three methods. Second, the impact of underground cable  #
# connections on the harmonic hosting capacity is investigated.                #
################################################################################
# Authors: Tom Van Acker, Hakan Ergun                                          #
################################################################################
# Changelog:                                                                   #
# v0.3.0 - init                                                                #
################################################################################

# INPUT ########################################################################
# using pkgs
using PowerModels
using HarmonicPowerModels

using Ipopt, JuMP

using Plots, StatsPlots, Measures

# pkg const
const PMs = PowerModels 
const HPM = HarmonicPowerModels

# solvers
solver_nlp = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)

# read in data
path = joinpath(HPM.BASE_DIR,"test/data/matpower/30bus_network_hhc.m")
data = PMs.parse_file(path)

# update the fairness principle
data["principle"] = "absolute equality"

# update the applicable standard of the buses
for (nb, bus) in data["bus"] 
    if bus["base_kv"] > 69.0
        bus["ref_angle"] = 0.0
        bus["standard"] = "IEEE519-2022-69/161kV"
    else
        bus["ref_angle"] = 0.0
        bus["standard"] = "IEEE519-2022-1/69kV"
    end
end

# update generator short-circuit impedances
Ssc     = [2.1e9, 2.1e9, 2.1e9, 2.1e9, 2.1e9, 2.1e9]
XRr     = 20
for (ng, gen) in data["gen"]
    gen["rsc"] = 100e6 / Ssc[parse(Int,ng)] / sqrt(1 + XRr^2)
    gen["xsc"] = 100e6 / Ssc[parse(Int,ng)] / sqrt(1 + (1/XRr)^2)
end

# update branch limits in base case to make case feasible
data["branch"]["22"]["rate_a"] = 0.25 * sqrt(3) / 0.9
data["branch"]["13"]["rate_a"] = 0.25 * sqrt(3) / 0.9

# update the data dictionary with i_load and i_ratio of the buses
## nodal three-phase short-circuit current [pu], determined using GDT.jl (private)
i_sc = [0.163, 14.94, 11.828, 12.944, 4.849, 10.76, 8.024, 5.19, 6.093, 6.09, 
        2.366, 5.342, 2.556, 3.55, 4.64, 4.234, 4.713, 3.426, 3.42, 3.669, 4.82, 
        4.806, 3.516, 4.033, 2.89, 1.312, 3.117, 8.267, 1.601, 1.461]
## init i_load and i_ratio for every bus
for (nb, bus) in data["bus"]
    bus["i_load"]  = 0.0 
    bus["i_ratio"] = Inf
end
## calculate i_load and i_ratio for all load buses
for (nl, load) in data["load"]
    nb = load["load_bus"]

    i_load = abs(load["pd"] + im * load["qd"])

    data["bus"]["$nb"]["i_load"]  = i_load
    data["bus"]["$nb"]["i_ratio"] = i_sc[nb] / i_load
end

# define the set of considered harmonics, and frequency
H = 1:50
f = 50:5:2500

# COMPUTATION ##################################################################
## A) Inductive Network ########################################################
# deepcopy of the original data
data_ind = deepcopy(data)

# reset all susceptance of the branches
for (nb, branch) in data_ind["branch"] 
    branch["b_fr"] = 0.0
    branch["b_to"] = 0.0
end

### IEEE519
# replicate the data
hdata_ind_519 = HPM.replicate(data_ind, H=collect(H))
# update the hdata with i_ihd and multiplier
for nh in setdiff(H,1)
    for (nb, bus) in hdata_ind_519["nw"]["$nh"]["bus"]
        v = bus["base_kv"]
        i = bus["i_ratio"]

        bus["i_ihd"] = HPM.ihd_current_limit_ieee_519(v, i, nh) * bus["i_load"] 
    end
    for (nl, load) in hdata_ind_519["nw"]["$nh"]["load"]
        nb = load["load_bus"]

        v = data["bus"]["$nb"]["base_kv"]
        i = data["bus"]["$nb"]["i_ratio"]

        load["multiplier"] = HPM.ihd_current_limit_ieee_519(v, i, nh)
end end
# solve the harmonic power flow
results_ind_519 = HPM.solve_hpf(hdata_ind_519, PMs.IVRPowerModel, solver_nlp)

### IEC61000-3-6
# replicate the data
hdata_ind_iec = HPM.replicate(data_ind, H=collect(H))
# calculate the impedance
Zh_ind = calculate_pos_seq_harmonic_impedance(data_ind, collect(f), collect(1:30))
# calculate the global contribution
gc = Dict(nb => bus["i_load"] / sum(bus["i_load"] for (nb,bus) in data["bus"]) for (nb,bus) in data["bus"])
# update the hdata with multiplier
for nh in setdiff(H,1)
    for (nb, bus) in hdata_ind_iec["nw"]["$nh"]["bus"]
        bus["i_ihd"] = 0.0
    end
    for (nl, load) in hdata_ind_iec["nw"]["$nh"]["load"]
        nb = load["load_bus"]
        ni = findfirst(x -> x == nh * 50, f)

        bus = hdata_ind_iec["nw"]["$nh"]["bus"]["$nb"]

        zh = abs(Zh_ind[nb][ni])
        uh = hdata_ind_iec["nw"]["$nh"]["bus"]["$nb"]["base_kv"] <= 69.0 ? 0.03 : 0.015 ;

        ihd = uh * gc["$nb"] / zh

        bus["i_ihd"] = ihd
        load["multiplier"] = ihd / bus["i_load"]
end end
# solve the harmonic power flow
results_ind_iec = HPM.solve_hpf(hdata_ind_iec, PMs.IVRPowerModel, solver_nlp)

### HHC Optimization Approach
hdata_ind_hhc   = HPM.replicate(data_ind, H=collect(H))
results_ind_hhc  = HPM.solve_hhc(hdata_ind_hhc, dHHC_SOC, solver_nlp, solver_nlp)

## B) Underground Cable Network ################################################
# deepcopy of the original data
data_cap = deepcopy(data)

### IEEE519
# replicate the data
hdata_cap_519 = HPM.replicate(data_cap, H=collect(H))
# update the hdata with i_ihd and multiplier
for nh in setdiff(H,1)
    for (nb, bus) in hdata_cap_519["nw"]["$nh"]["bus"]
        v = bus["base_kv"]
        i = bus["i_ratio"]

        bus["i_ihd"] = HPM.ihd_current_limit_ieee_519(v, i, nh) * bus["i_load"] 
    end
    for (nl, load) in hdata_cap_519["nw"]["$nh"]["load"]
        nb = load["load_bus"]

        v = data["bus"]["$nb"]["base_kv"]
        i = data["bus"]["$nb"]["i_ratio"]

        load["multiplier"] = HPM.ihd_current_limit_ieee_519(v, i, nh)
end end
# solve the harmonic power flow
results_cap_519 = HPM.solve_hpf(hdata_cap_519, PMs.IVRPowerModel, solver_nlp)

### IEC61000-3-6
# replicate the data
hdata_cap_iec = HPM.replicate(data_cap, H=collect(H))
# calculate the impedance
Zh_cap = calculate_pos_seq_harmonic_impedance(data_cap, collect(f), collect(1:30))
# calculate the global contribution
gc = Dict(nb => bus["i_load"] / sum(bus["i_load"] for (nb,bus) in data["bus"]) for (nb,bus) in data["bus"])
# update the hdata with multiplier
for nh in setdiff(H,1)
    for (nb, bus) in hdata_cap_iec["nw"]["$nh"]["bus"]
        bus["i_ihd"] = 0.0
    end
    for (nl, load) in hdata_cap_iec["nw"]["$nh"]["load"]
        nb = load["load_bus"]
        ni = findfirst(x -> x == nh * 50, f)

        bus = hdata_cap_iec["nw"]["$nh"]["bus"]["$nb"]

        zh = abs(Zh_cap[nb][ni])
        uh = hdata_cap_iec["nw"]["$nh"]["bus"]["$nb"]["base_kv"] <= 69.0 ? 0.03 : 0.015 ;

        ihd = uh * gc["$nb"] / zh

        bus["i_ihd"] = ihd
        load["multiplier"] = ihd / bus["i_load"]
end end
# solve the harmonic power flow
results_cap_iec = HPM.solve_hpf(hdata_cap_iec, PMs.IVRPowerModel, solver_nlp)

### HHC Optimization Approach
hdata_cap_hhc   = HPM.replicate(data_cap, H=collect(H))
results_cap_hhc  = HPM.solve_hhc(hdata_cap_hhc, dHHC_SOC, solver_nlp, solver_nlp)

# RESULTS ######################################################################
# select a load for which to plot results
nl = 3
nb = data["load"]["$nl"]["load_bus"]

# plot harmonic impedance
plot(f ./ 50, [abs.(Zh_ind[nb]), abs.(Zh_cap[nb])],
        fontfamily="Computer Modern",
        label=["Case A: Inductive network" "Case B: Underground cable network"],
        legend=:bottomright,
        linestyle=[:solid :dash],
        xlabel="harmonic \$h\$ [-]",
        yaxis=:log10,
        ylabel="impedance \$Z_{h}\$ [pu]",
        ylim=(1e-2,1e1),
        yminorgrid=true)
savefig("examples/dhhc_pes_gm/results/harmonic_impedance_bus_$nb.png")

# plot harmonic current
plot(setdiff(H,1), [[hdata_ind_519["nw"]["$nh"]["bus"]["$nb"]["i_ihd"] for nh in setdiff(H,1)],
                    [hdata_ind_iec["nw"]["$nh"]["bus"]["$nb"]["i_ihd"] for nh in setdiff(H,1)],
                    [results_ind_hhc["solution"]["nw"]["$nh"]["load"]["$nl"]["cmd"] for nh in setdiff(H,1)],
                    [hdata_cap_iec["nw"]["$nh"]["bus"]["$nb"]["i_ihd"] for nh in setdiff(H,1)],
                    [results_cap_hhc["solution"]["nw"]["$nh"]["load"]["$nl"]["cmd"] for nh in setdiff(H,1)]
                   ],
        bottom_margin=10mm,
        color=[:blue :orange :green :orange :green],
        fontfamily="Computer Modern",
        label=["CASE A/B: IEEE519" "CASE B: IEC61000" "CASE B: HHC" "CASE A: IEC61000" "CASE A: HHC"],
        left_margin=10mm,
        legend=:outerright,
        linestyle=[:dot :dot :dot :solid :solid :solid],  
        marker=[:square :square :square :diamond :diamond :diamond],  
        markersize=2,
        size=(1200,266),
        xlabel="harmonic \$h\$ [-]",
        xlim=(0,51),
        yaxis=:log10,
        ylabel="current \$|𝗜^{ihd}_{h}|\$ [pu]",
        ylim=(1e-5,1e-2),
        yminorgrid=true)
savefig("examples/dhhc_pes_gm/results/harmonic_current_load_$nl.png")

nb = 26

# plot harmonic voltage
vm(result) = abs(result["vr"] + im * result["vi"])
lim = data["bus"]["$nb"]["base_kv"] <= 69.0 ? 0.03 : 0.015 ;
plot([0,50],[lim,lim],color=:black,label="harmonic voltage limit",linewidth=3)
plot!(setdiff(H,1), [[vm(results_ind_519["solution"]["nw"]["$nh"]["bus"]["$nb"]) for nh in setdiff(H,1)],
                     [vm(results_ind_iec["solution"]["nw"]["$nh"]["bus"]["$nb"]) for nh in setdiff(H,1)],
                     [vm(results_ind_hhc["solution"]["nw"]["$nh"]["bus"]["$nb"]) for nh in setdiff(H,1)],
                     [vm(results_cap_519["solution"]["nw"]["$nh"]["bus"]["$nb"]) for nh in setdiff(H,1)],
                     [vm(results_cap_iec["solution"]["nw"]["$nh"]["bus"]["$nb"]) for nh in setdiff(H,1)],
                     [vm(results_cap_hhc["solution"]["nw"]["$nh"]["bus"]["$nb"]) for nh in setdiff(H,1)]
                    ],
        bottom_margin=10mm,
        color=[:blue :orange :green :blue :orange :green],
        fontfamily="Computer Modern",
        label=["CASE A: IEEE519" "CASE A: IEC61000" "CASE A: HHC" "CASE B: IEEE519" "CASE B: IEC61000" "CASE B: HHC"],  
        left_margin=10mm,
        legend=:outerright,
        linestyle=[:solid :solid :solid :dash :dash :dash], 
        marker=[:diamond :diamond :diamond :square :square :square],  
        markersize=5,
        size=(1200,266),
        xlabel="harmonic \$h\$ [-]",
        xlim=(0,51),
        yaxis=:log10,
        ylabel="voltage \$|𝗨^{ihd}_{h}|\$ [pu]",
        ylim=(1e-5,1e0),
        yminorgrid=true)
savefig("examples/dhhc_pes_gm/results/harmonic_voltage_bus_$nb.png")


######### IHD IND

nb = 4
# plot harmonic current
 p_indhv=plot(setdiff(H,1), [[hdata_ind_519["nw"]["$nh"]["bus"]["$nb"]["i_ihd"] for nh in setdiff(H,1)],
                    [hdata_ind_iec["nw"]["$nh"]["bus"]["$nb"]["i_ihd"] for nh in setdiff(H,1)],
                    [results_ind_hhc["solution"]["nw"]["$nh"]["load"]["$nl"]["cmd"] for nh in setdiff(H,1)]
                    # [hdata_cap_iec["nw"]["$nh"]["bus"]["$nb"]["i_ihd"] for nh in setdiff(H,1)],
                    # [results_cap_hhc["solution"]["nw"]["$nh"]["load"]["$nl"]["cmd"] for nh in setdiff(H,1)]
                   ],
        bottom_margin=10mm,
        color=[:blue :orange :green],
        fontfamily="Computer Modern",
        label=["IEEE519" "IEC61000" "HHC"],
        left_margin=10mm,
        legend=:outerright,
        linestyle=[:dot :dot :dot],  
        marker=[:square :square :square],  
        markersize=2,
        size=(1200,266),
        xlabel="harmonic \$h\$ [-]",
        xlim=(0,51),
        yaxis=:log10,
        ylabel="current \$|𝗜^{ihd}_{h}|\$ [pu]",
        ylim=(1e-5,1e-1),
        guidefontsize=16,
        tickfontsize=14,
        legendfontsize=16,
        linewidth = 3,
        yminorgrid=true)
savefig(p_indhv, "ihd_ind_hv.pdf")

nb = 17
# plot harmonic current
p_indmv= plot(setdiff(H,1), [[hdata_ind_519["nw"]["$nh"]["bus"]["$nb"]["i_ihd"] for nh in setdiff(H,1)],
                    [hdata_ind_iec["nw"]["$nh"]["bus"]["$nb"]["i_ihd"] for nh in setdiff(H,1)],
                    [results_ind_hhc["solution"]["nw"]["$nh"]["load"]["$nl"]["cmd"] for nh in setdiff(H,1)]
                    # [hdata_cap_iec["nw"]["$nh"]["bus"]["$nb"]["i_ihd"] for nh in setdiff(H,1)],
                    # [results_cap_hhc["solution"]["nw"]["$nh"]["load"]["$nl"]["cmd"] for nh in setdiff(H,1)]
                   ],
         bottom_margin=10mm,
        color=[:blue :orange :green],
        fontfamily="Computer Modern",
        label=["IEEE519" "IEC61000" "HHC"],
        left_margin=10mm,
        legend=:outerright,
        linestyle=[:dot :dot :dot],  
        marker=[:square :square :square],  
        markersize=2,
        size=(1200,266),
        xlabel="harmonic \$h\$ [-]",
        xlim=(0,51),
        yaxis=:log10,
        ylabel="current \$|𝗜^{ihd}_{h}|\$ [pu]",
        ylim=(1e-5,1e-1),
        guidefontsize=16,
        tickfontsize=14,
        legendfontsize=16,
        linewidth = 3,
        yminorgrid=true)
savefig(p_indmv, "ihd_ind_mv.pdf")


# plot harmonic voltage
nb = 4
vm(result) = abs(result["vr"] + im * result["vi"])
lim = data["bus"]["$nb"]["base_kv"] <= 69.0 ? 0.03 : 0.015 ;
pv_ind_hv = plot([0,50],[lim,lim],color=:black,label="harmonic voltage limit",linewidth=3)
plot!(setdiff(H,1), [[vm(results_ind_519["solution"]["nw"]["$nh"]["bus"]["$nb"]) for nh in setdiff(H,1)],
                     [vm(results_ind_iec["solution"]["nw"]["$nh"]["bus"]["$nb"]) for nh in setdiff(H,1)],
                     [vm(results_ind_hhc["solution"]["nw"]["$nh"]["bus"]["$nb"]) for nh in setdiff(H,1)]
                    #  [vm(results_cap_519["solution"]["nw"]["$nh"]["bus"]["$nb"]) for nh in setdiff(H,1)],
                    #  [vm(results_cap_iec["solution"]["nw"]["$nh"]["bus"]["$nb"]) for nh in setdiff(H,1)],
                    #  [vm(results_cap_hhc["solution"]["nw"]["$nh"]["bus"]["$nb"]) for nh in setdiff(H,1)]
                    ],
         bottom_margin=10mm,
        color=[:blue :orange :green],
        fontfamily="Computer Modern",
        label=["IEEE519" "IEC61000" "HHC"],
        left_margin=10mm,
        legend=:outerright,
        linestyle=[:dot :dot :dot],  
        marker=[:square :square :square],  
        markersize=2,
        size=(1200,266),
        xlabel="harmonic \$h\$ [-]",
        xlim=(0,51),
        yaxis=:log10,
        ylabel="voltage \$|𝗨^{ihd}_{h}|\$ [pu]",
        ylim=(1e-5,1e-1),
        guidefontsize=16,
        tickfontsize=14,
        legendfontsize=16,
        linewidth = 3,
        yminorgrid=true)
savefig(pv_ind_hv, "v_ind_hv.pdf")


# plot harmonic voltage
nb = 17
vm(result) = abs(result["vr"] + im * result["vi"])
lim = data["bus"]["$nb"]["base_kv"] <= 69.0 ? 0.03 : 0.015 ;
pv_ind_mv = plot([0,50],[lim,lim],color=:black,label="harmonic voltage limit",linewidth=3)
plot!(setdiff(H,1), [[vm(results_ind_519["solution"]["nw"]["$nh"]["bus"]["$nb"]) for nh in setdiff(H,1)],
                     [vm(results_ind_iec["solution"]["nw"]["$nh"]["bus"]["$nb"]) for nh in setdiff(H,1)],
                     [vm(results_ind_hhc["solution"]["nw"]["$nh"]["bus"]["$nb"]) for nh in setdiff(H,1)]
                    #  [vm(results_cap_519["solution"]["nw"]["$nh"]["bus"]["$nb"]) for nh in setdiff(H,1)],
                    #  [vm(results_cap_iec["solution"]["nw"]["$nh"]["bus"]["$nb"]) for nh in setdiff(H,1)],
                    #  [vm(results_cap_hhc["solution"]["nw"]["$nh"]["bus"]["$nb"]) for nh in setdiff(H,1)]
                    ],
         bottom_margin=10mm,
        color=[:blue :orange :green],
        fontfamily="Computer Modern",
        label=["IEEE519" "IEC61000" "HHC"],
        left_margin=10mm,
        legend=:outerright,
        linestyle=[:dot :dot :dot],  
        marker=[:square :square :square],  
        markersize=2,
        size=(1200,266),
        xlabel="harmonic \$h\$ [-]",
        xlim=(0,51),
        yaxis=:log10,
        ylabel="voltage \$|𝗨^{ihd}_{h}|\$ [pu]",
        ylim=(1e-5,1e-1),
        guidefontsize=16,
        tickfontsize=14,
        legendfontsize=16,
        linewidth = 3,
        yminorgrid=true)
savefig(pv_ind_mv, "v_ind_mv.pdf")




######### IHD CAP

nb = 4
# plot harmonic current
 p_caphv = plot(setdiff(H,1), [[hdata_ind_519["nw"]["$nh"]["bus"]["$nb"]["i_ihd"] for nh in setdiff(H,1)],
                    [hdata_cap_iec["nw"]["$nh"]["bus"]["$nb"]["i_ihd"] for nh in setdiff(H,1)],
                    [results_cap_hhc["solution"]["nw"]["$nh"]["load"]["$nl"]["cmd"] for nh in setdiff(H,1)]
                    # [hdata_cap_iec["nw"]["$nh"]["bus"]["$nb"]["i_ihd"] for nh in setdiff(H,1)],
                    # [results_cap_hhc["solution"]["nw"]["$nh"]["load"]["$nl"]["cmd"] for nh in setdiff(H,1)]
                   ],
        bottom_margin=10mm,
        color=[:blue :orange :green],
        fontfamily="Computer Modern",
        label=["IEEE519" "IEC61000" "HHC"],
        left_margin=10mm,
        legend=:outerright,
        linestyle=[:dot :dot :dot],  
        marker=[:square :square :square],  
        markersize=2,
        size=(1200,266),
        xlabel="harmonic \$h\$ [-]",
        xlim=(0,51),
        yaxis=:log10,
        ylabel="current \$|𝗜^{ihd}_{h}|\$ [pu]",
        ylim=(1e-5,1e-1),
        guidefontsize=16,
        tickfontsize=14,
        legendfontsize=16,
        linewidth = 3,
        yminorgrid=true)
savefig(p_caphv, "ihd_cap_hv.pdf")

nb = 17
# plot harmonic current
p_capmv= plot(setdiff(H,1), [[hdata_ind_519["nw"]["$nh"]["bus"]["$nb"]["i_ihd"] for nh in setdiff(H,1)],
                    [hdata_cap_iec["nw"]["$nh"]["bus"]["$nb"]["i_ihd"] for nh in setdiff(H,1)],
                    [results_cap_hhc["solution"]["nw"]["$nh"]["load"]["$nl"]["cmd"] for nh in setdiff(H,1)]
                    # [hdata_cap_iec["nw"]["$nh"]["bus"]["$nb"]["i_ihd"] for nh in setdiff(H,1)],
                    # [results_cap_hhc["solution"]["nw"]["$nh"]["load"]["$nl"]["cmd"] for nh in setdiff(H,1)]
                   ],
         bottom_margin=10mm,
        color=[:blue :orange :green],
        fontfamily="Computer Modern",
        label=["IEEE519" "IEC61000" "HHC"],
        left_margin=10mm,
        legend=:outerright,
        linestyle=[:dot :dot :dot],  
        marker=[:square :square :square],  
        markersize=2,
        size=(1200,266),
        xlabel="harmonic \$h\$ [-]",
        xlim=(0,51),
        yaxis=:log10,
        ylabel="current \$|𝗜^{ihd}_{h}|\$ [pu]",
        ylim=(1e-5,1e-1),
        guidefontsize=16,
        tickfontsize=14,
        legendfontsize=16,
        linewidth = 3,
        yminorgrid=true)
savefig(p_capmv, "ihd_cap_mv.pdf")



# plot harmonic voltage
nb = 4
vm(result) = abs(result["vr"] + im * result["vi"])
lim = data["bus"]["$nb"]["base_kv"] <= 69.0 ? 0.03 : 0.015 ;
pv_cap_hv = plot([0,50],[lim,lim],color=:black,label="harmonic voltage limit",linewidth=3)
plot!(setdiff(H,1), [
                     [vm(results_cap_519["solution"]["nw"]["$nh"]["bus"]["$nb"]) for nh in setdiff(H,1)],
                     [vm(results_cap_iec["solution"]["nw"]["$nh"]["bus"]["$nb"]) for nh in setdiff(H,1)],
                     [vm(results_cap_hhc["solution"]["nw"]["$nh"]["bus"]["$nb"]) for nh in setdiff(H,1)]
                    ],
         bottom_margin=10mm,
        color=[:blue :orange :green],
        fontfamily="Computer Modern",
        label=["IEEE519" "IEC61000" "HHC"],
        left_margin=10mm,
        legend=:outerright,
        linestyle=[:dot :dot :dot],  
        marker=[:square :square :square],  
        markersize=2,
        size=(1200,266),
        xlabel="harmonic \$h\$ [-]",
        xlim=(0,51),
        yaxis=:log10,
        ylabel="voltage \$|𝗨^{ihd}_{h}|\$ [pu]",
        ylim=(1e-5,1e-0),
        guidefontsize=16,
        tickfontsize=14,
        legendfontsize=16,
        linewidth = 3,
        yminorgrid=true)
savefig(pv_cap_hv, "v_cap_hv.pdf")


# plot harmonic voltage
nb = 17
vm(result) = abs(result["vr"] + im * result["vi"])
lim = data["bus"]["$nb"]["base_kv"] <= 69.0 ? 0.03 : 0.015 ;
pv_cap_mv = plot([0,50],[lim,lim],color=:black,label="harmonic voltage limit",linewidth=3)
plot!(setdiff(H,1), [
                     [vm(results_cap_519["solution"]["nw"]["$nh"]["bus"]["$nb"]) for nh in setdiff(H,1)],
                     [vm(results_cap_iec["solution"]["nw"]["$nh"]["bus"]["$nb"]) for nh in setdiff(H,1)],
                     [vm(results_cap_hhc["solution"]["nw"]["$nh"]["bus"]["$nb"]) for nh in setdiff(H,1)]
                    ],
         bottom_margin=10mm,
        color=[:blue :orange :green],
        fontfamily="Computer Modern",
        label=["IEEE519" "IEC61000" "HHC"],
        left_margin=10mm,
        legend=:outerright,
        linestyle=[:dot :dot :dot],  
        marker=[:square :square :square],  
        markersize=2,
        size=(1200,266),
        xlabel="harmonic \$h\$ [-]",
        xlim=(0,51),
        yaxis=:log10,
        ylabel="voltage \$|𝗨^{ihd}_{h}|\$ [pu]",
        ylim=(1e-5,1e-0),
        guidefontsize=16,
        tickfontsize=14,
        legendfontsize=16,
        linewidth = 3,
        yminorgrid=true)
savefig(pv_cap_mv, "v_cap_mv.pdf")


nb = 4
nb_1 = 17
# plot harmonic impedance
p1 = plot(f ./ 50, [abs.(Zh_ind[nb])],
        fontfamily="Computer Modern",
        label="HV node",
        legend=:bottomright,
        color = RGB(0,102/255,51/255),
        linestyle=:solid,
        xlabel="harmonic \$h\$ [-]",
        yaxis=:log10,
        ylabel="impedance \$Z_{h}\$ [pu]",
        ylim=(1e-2,1e1),
        linewidth = 4,
        guidefontsize=16,
        tickfontsize=14,
        legendfontsize=16,
        yminorgrid=true)




plot!(p1, f ./ 50, [abs.(Zh_ind[nb_1])],
    fontfamily="Computer Modern",
    label="MV node",
    legend=:bottomright,
        color = RGB(0,102/255,51/255),
    linestyle=:dot,
    xlabel="harmonic \$h\$ [-]",
    yaxis=:log10,
    ylabel="impedance \$Z_{h}\$ [pu]",
    ylim=(1e-2,1e1),
    linewidth = 4,
    guidefontsize=16,
    tickfontsize=14,
    legendfontsize=16,
    yminorgrid=true)


savefig(p1, "zh_ind.pdf")



nb = 4
nb_1 = 17
# plot harmonic impedance
p2 = plot(f ./ 50, [abs.(Zh_cap[nb])],
        fontfamily="Computer Modern",
        label="HV node",
        legend=:bottomright,
        color = RGB(0,102/255,51/255),
        linestyle=:solid,
        xlabel="harmonic \$h\$ [-]",
        yaxis=:log10,
        ylabel="impedance \$Z_{h}\$ [pu]",
        ylim=(1e-2,1e1),
        linewidth = 4,
        guidefontsize=16,
        tickfontsize=14,
        legendfontsize=16,
        yminorgrid=true)


plot!(p2, f ./ 50, [abs.(Zh_cap[nb_1])],
    fontfamily="Computer Modern",
    label="MV node",
    legend=:bottomright,
        color = RGB(0,102/255,51/255),
    linestyle=:dot,
    xlabel="harmonic \$h\$ [-]",
    yaxis=:log10,
    ylabel="impedance \$Z_{h}\$ [pu]",
    ylim=(1e-2,1e1),
    linewidth = 4,
    guidefontsize=16,
    tickfontsize=14,
    legendfontsize=16,
    yminorgrid=true)


savefig(p2, "zh_cap.pdf")