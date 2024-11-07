# using pkgs
using PowerModels
using HarmonicPowerModels

using Ipopt, JuMP

using Plots

using Revise

# pkg const
const PMs = PowerModels 
const HPM = HarmonicPowerModels

# read in data
path = joinpath(HPM.BASE_DIR,"test/data/matpower/30bus_network_hhc_519.m")
data = PMs.parse_file(path)

# solvers
solver_nlp = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)

# reset all susceptance of the branches
for (nb, branch) in data["branch"] 
    branch["b_fr"] = 0.0
    branch["b_to"] = 0.0
end

# add the load current magnitude [pu] to all buses
for (nb, bus) in data["bus"]
    bus["i_load"] = 0.0
end
for (nl, load) in data["load"]
    nb = load["load_bus"]
    data["bus"]["$nb"]["i_load"] = abs(load["pd"] + im * load["qd"])
end

# add the short-circuit current magnitude [pu] to all busses
i_sc = [0.163, 14.94, 11.828, 12.944, 4.849, 10.76, 8.024, 5.19, 6.093, 6.09, 
        2.366, 5.342, 2.556, 3.55, 4.64, 4.234, 4.713, 3.426, 3.42, 3.669, 4.82, 
        4.806, 3.516, 4.033, 2.89, 1.312, 3.117, 8.267, 1.601, 1.461]
for (nb, bus) in data["bus"]
    bus["i_sc"] = i_sc[parse(Int, nb)]
end

# add the short-circuit vs load current ratio
for (nb, bus) in data["bus"]
    bus["i_ratio"] = bus["i_sc"] / bus["i_load"]
end

# III.A purely inductive network ###############################################

# IEEE519-2022
## determine the individual harmonic unit currents
function ihd_limit_519(voltage, i_ratio, harmonic)
    if voltage < 69.0
        if i_ratio == Inf
            return 0.0
        elseif i_ratio < 20
            if 2 <= harmonic < 11
                return 0.04
            elseif 11 <= harmonic < 17
                return 0.02
            elseif 17 <= harmonic < 23
                return 0.015
            elseif 23 <= harmonic < 35
                return 0.006
            elseif 35 <= harmonic <= 50
                return 0.003
            else
                return 0.0
            end
        elseif 20.0 <= i_ratio < 50.0
            if 2 <= harmonic < 11
                return 0.07
            elseif 11 <= harmonic < 17
                return 0.035
            elseif 17 <= harmonic < 23
                return 0.025
            elseif 23 <= harmonic < 35
                return 0.010
            elseif 35 <= harmonic <= 50
                return 0.005
            else
                return 0.0
            end
        elseif 50.0 <= i_ratio < 100.0
            if 2 <= harmonic < 11
                return 0.10
            elseif 11 <= harmonic < 17
                return 0.045
            elseif 17 <= harmonic < 23
                return 0.040
            elseif 23 <= harmonic < 35
                return 0.015
            elseif 35 <= harmonic <= 50
                return 0.007
            else
                return 0.0
            end
        elseif 100.0 <= i_ratio < 1000.0
            if 2 <= harmonic < 11
                return 0.12
            elseif 11 <= harmonic < 17
                return 0.055
            elseif 17 <= harmonic < 23
                return 0.050
            elseif 23 <= harmonic < 35
                return 0.020
            elseif 35 <= harmonic <= 50
                return 0.010
            else
                return 0.0
            end
        else
            if 2 <= harmonic < 11
                return 0.15
            elseif 11 <= harmonic < 17
                return 0.07
            elseif 17 <= harmonic < 23
                return 0.06
            elseif 23 <= harmonic < 35
                return 0.025
            elseif 35 <= harmonic <= 50
                return 0.014
            else
                return 0.0
            end
        end
    else
        if i_ratio == Inf
            return 0.0
        elseif i_ratio < 20       
            if 2 <= harmonic < 11
                return 0.02
            elseif 11 <= harmonic < 17
                return 0.01
            elseif 17 <= harmonic < 23
                return 0.0075
            elseif 23 <= harmonic < 35
                return 0.003
            elseif 35 <= harmonic <= 50
                return 0.0015
            else
                return 0.0
            end
        elseif 20.0 <= i_ratio < 50.0
            if 2 <= harmonic < 11
                return 0.035
            elseif 11 <= harmonic < 17
                return 0.0175
            elseif 17 <= harmonic < 23
                return 0.0125
            elseif 23 <= harmonic < 35
                return 0.005
            elseif 35 <= harmonic <= 50
                return 0.0025
            else
                return 0.0
            end
        elseif 50.0 <= i_ratio < 100.0
            if 2 <= harmonic < 11
                return 0.05
            elseif 11 <= harmonic < 17
                return 0.0225
            elseif 17 <= harmonic < 23
                return 0.02
            elseif 23 <= harmonic < 35
                return 0.0075
            elseif 35 <= harmonic <= 50
                return 0.0035
            else
                return 0.0
            end
        elseif 100.0 <= i_ratio < 1000.0
            if 2 <= harmonic < 11
                return 0.06
            elseif 11 <= harmonic < 17
                return 0.0275
            elseif 17 <= harmonic < 23
                return 0.025
            elseif 23 <= harmonic < 35
                return 0.01
            elseif 35 <= harmonic <= 50
                return 0.005
            else
                return 0.0
            end
        else
            if 2 <= harmonic < 11
                return 0.075
            elseif 11 <= harmonic < 17
                return 0.035
            elseif 17 <= harmonic < 23
                return 0.03
            elseif 23 <= harmonic < 35
                return 0.0125
            elseif 35 <= harmonic <= 50
                return 0.007
            else
                return 0.0
            end
        end
    end
end
for (nb, bus) in data["bus"]
    v = bus["base_kv"]
    i = bus["i_ratio"]
    bus["i_ihd_519"] = [ihd_limit_519(v, i, h) * bus["i_load"] for h in 2:50]
end

## solve harmonic load flow to determine the bus voltage
hdata_ind_519       = HPM.replicate(data, H=collect(1:50))
for (nw, ntw) in hdata_ind_519["nw"] if parse(Int, nw) ≠ 1
    for (nl, load) in ntw["load"]
        nb = load["load_bus"]

        v = data["bus"]["$nb"]["base_kv"]
        i = data["bus"]["$nb"]["i_ratio"]

        load["multiplier"] = ihd_limit_519(v, i , parse(Int,nw))
end end end
results_hpf_ind_519 = HPM.solve_hpf(hdata_ind_519, PMs.IVRPowerModel, solver_nlp)


# Plots
## input
nb = 17
nl = 1
## function
vm(bus)     = abs(bus["vr"] + im * bus["vi"])

## plots 
bar(2:50, [data["bus"]["$nb"]["i_ihd_519"][nh-1] for nh in 2:50], 
        label="IEEE519-2022",
        xlabel="harmonic h [-]",
        ylabel="individual harmonic unit current Iⁱʰᵈ [pu]")
bar(2:50, [results_hpf_ind_519["solution"]["nw"]["$nh"]["load"]["$nl"]["cm"] for nh in 2:50])
bar(1:50, [vm(results_hpf_ind_519["solution"]["nw"]["$nh"]["bus"]["$nb"]) for nh in 1:50],
        label="IEEE519-2022",
        xlabel="harmonic h [-]",
        ylabel="individual harmonic bus voltage Uⁱʰᵈ [pu]")