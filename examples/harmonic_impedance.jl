#######  Analysis for PES GM paper  #####

using HarmonicPowerModels
using PowerModels
using Plots
using Ipopt
using StatsPlots
using Measures

const PMs = PowerModels
const HPM = HarmonicPowerModels

## Include function to calculate Harmonic impedance
include(joinpath(HPM.BASE_DIR,"src/util/hi.jl"))

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

# Define function to plot voltages
vm(bus)     = abs(bus["vr"] + im * bus["vi"])

## Load test case
path = joinpath(HPM.BASE_DIR,"test/data/matpower/30bus_network_hhc_519.m")
data = PMs.parse_file(path)


# set the solver
solver_nlp = Ipopt.Optimizer

# reset all susceptance of the branches
for (nb, branch) in data["branch"] 
    branch["b_fr"] = 0.0
    branch["b_to"] = 0.0
end

# reset all susceptance of the branches
for (nb, bus) in data["bus"] 
    if bus["base_kv"] > 69.0
        bus["ref_angle"] = 0.0
        bus["standard"] = "IEEE519-2022-69/161kV"
    else
        bus["ref_angle"] = 0.0
        bus["standard"] = "IEEE519-2022-1/69kV"
    end
end

for (nl, load) in data["load"]
    load["i_load"] = abs(load["pd"] + im * load["qd"])
end

# change branch limits in base case to make case feasible
data["branch"]["22"]["rate_a"] = 0.25 * sqrt(3) / 0.9
data["branch"]["13"]["rate_a"] = 0.25 * sqrt(3) / 0.9


###################################################
###### Determine hosting capacity with optimisation

# define the set of considered harmonics
H = [i for i in 1:50]

# Principle -> Absolute equality
data["principle"] = "absolute equality"

# solve HHC problem -- SOC
hdata = HPM.replicate(data, H=H)

# Optimisation
result_hhc = HPM.solve_hhc(hdata, dHHC_SOC, solver_nlp, solver_nlp)

###################################################
###### Determine hosting with IEC61000

# Step 1: Calculate Harmonic impedance
Ssc     = 2.1e9
Sbase   = 100e6
XRr     = 20.0

z = Sbase / Ssc
r = z / sqrt(1 + XRr^2)
x = sqrt(z^2 - r^2)
for (ng, gen) in data["gen"]
    gen["rsc"] = r
    gen["xsc"] = x
end

# define parameters for harmonic impedance calculation
number_of_buses = 30
Hf = 50:5.0:2500

# Calculate impedance
Zh = calculate_pos_seq_harmonic_impedance(data, collect(Hf), collect(1:number_of_buses))

# Calculate the global contribution per load bus
global_contribution = [abs.(data["load"]["$l"]["pd"] + im * data["load"]["$l"]["qd"]) for l in sort(parse.(Int, keys(data["load"])))] ./ sum([abs.(load["pd"] + im * load["qd"]) for (l, load) in data["load"]])

hdata = HPM.replicate(data, H=H)
for (nw, ntw) in hdata["nw"]
    for (l, load) in ntw["load"]
        bus_id = load["load_bus"]
        h_id = (parse(Int, nw) - 1) * 10 + 1
        zh = Zh[bus_id][h_id]
        if ntw["bus"]["$bus_id"]["base_kv"] <= 69.0
            ihd = 0.03 * global_contribution[parse(Int, l)] / abs(zh)
        else
            ihd = 0.015 * global_contribution[parse(Int, l)] / abs(zh)
        end
        load["multiplier"] = ihd / load["i_load"]
    end
end

# Solve Harmonic power flow
result_iec = HPM.solve_hpf(hdata, PMs.IVRPowerModel, solver_nlp)


#######################################################
########## Calcukate harmonic currents based on IEEE 519
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
result_ieee = HPM.solve_hpf(hdata_ind_519, PMs.IVRPowerModel, solver_nlp)

################################################################
###########     UNDERGROUND CABLE   ###################################
################################################################
################################################################
path = joinpath(HPM.BASE_DIR,"test/data/matpower/30bus_network_hhc_519.m")
data_cable = PMs.parse_file(path)
# reset all susceptance of the branches
for (nb, bus) in data_cable["bus"] 
    if bus["base_kv"] > 69.0
        bus["ref_angle"] = 0.0
        bus["standard"] = "IEEE519-2022-69/161kV"
    else
        bus["ref_angle"] = 0.0
        bus["standard"] = "IEEE519-2022-1/69kV"
    end
end

for (nl, load) in data_cable["load"]
    load["i_load"] = abs(load["pd"] + im * load["qd"])
end

# change branch limits in base case to make case feasible
data_cable["branch"]["22"]["rate_a"] = 0.25 * sqrt(3) / 0.9
data_cable["branch"]["13"]["rate_a"] = 0.25 * sqrt(3) / 0.9


###################################################
###### Determine hosting capacity with optimisation

# define the set of considered harmonics
H = [i for i in 1:50]

# Principle -> Absolute equality
data_cable["principle"] = "absolute equality"

# solve HHC problem -- SOC
hdata_cable = HPM.replicate(data_cable, H=H)

# Optimisation
result_hhc_cable = HPM.solve_hhc(hdata_cable, dHHC_SOC, solver_nlp, solver_nlp)

###################################################
###### Determine hosting with IEC61000

# Step 1: Calculate Harmonic impedance
Ssc     = 2.1e9
Sbase   = 100e6
XRr     = 20.0

z = Sbase / Ssc
r = z / sqrt(1 + XRr^2)
x = sqrt(z^2 - r^2)
for (ng, gen) in data_cable["gen"]
    gen["rsc"] = r
    gen["xsc"] = x
end

# define parameters for harmonic impedance calculation
number_of_buses = 30
Hf = 50:5.0:2500

# Calculate impedance
Zh_cable = calculate_pos_seq_harmonic_impedance(data_cable, collect(Hf), collect(1:number_of_buses))

# Calculate the global contribution per load bus
global_contribution = [abs.(data_cable["load"]["$l"]["pd"] + im * data_cable["load"]["$l"]["qd"]) for l in sort(parse.(Int, keys(data_cable["load"])))] ./ sum([abs.(load["pd"] + im * load["qd"]) for (l, load) in data_cable["load"]])

# Calcuate harmonic current


hdata_cable = HPM.replicate(data_cable, H=H)
for (nw, ntw) in hdata_cable["nw"]
    for (l, load) in ntw["load"]
        bus_id = load["load_bus"]
        h_id = (parse(Int, nw) - 1) * 10 + 1
        zh = Zh_cable[bus_id][h_id]
        if ntw["bus"]["$bus_id"]["base_kv"] <= 69.0
            ihd = 0.03 * global_contribution[parse(Int, l)] / abs(zh)
        else
            ihd = 0.015 * global_contribution[parse(Int, l)] / abs(zh)
        end
        load["multiplier"] = ihd / load["i_load"]
    end
end

# Solve Harmonic power flow
result_iec_cable = HPM.solve_hpf(hdata_cable, PMs.IVRPowerModel, solver_nlp)


#######################################################
########## Calcukate harmonic currents based on IEEE 519
# add the load current magnitude [pu] to all buses
for (nb, bus) in data_cable["bus"]
    bus["i_load"] = 0.0
end
for (nl, load) in data_cable["load"]
    nb = load["load_bus"]
    data_cable["bus"]["$nb"]["i_load"] = abs(load["pd"] + im * load["qd"])
end

# add the short-circuit current magnitude [pu] to all busses
i_sc = [0.163, 14.94, 11.828, 12.944, 4.849, 10.76, 8.024, 5.19, 6.093, 6.09, 
        2.366, 5.342, 2.556, 3.55, 4.64, 4.234, 4.713, 3.426, 3.42, 3.669, 4.82, 
        4.806, 3.516, 4.033, 2.89, 1.312, 3.117, 8.267, 1.601, 1.461]
for (nb, bus) in data_cable["bus"]
    bus["i_sc"] = i_sc[parse(Int, nb)]
end

# add the short-circuit vs load current ratio
for (nb, bus) in data_cable["bus"]
    bus["i_ratio"] = bus["i_sc"] / bus["i_load"]
end

# III.A purely inductive network ###############################################


for (nb, bus) in data_cable["bus"]
    v = bus["base_kv"]
    i = bus["i_ratio"]
    bus["i_ihd_519"] = [ihd_limit_519(v, i, h) * bus["i_load"] for h in 2:50]
end

## solve harmonic load flow to determine the bus voltage
hdata_ieee_cable = HPM.replicate(data_cable, H=collect(1:50))
for (nw, ntw) in hdata_ieee_cable["nw"] if parse(Int, nw) ≠ 1
    for (nl, load) in ntw["load"]
        nb = load["load_bus"]

        v = data["bus"]["$nb"]["base_kv"]
        i = data["bus"]["$nb"]["i_ratio"]

        load["multiplier"] = ihd_limit_519(v, i , parse(Int,nw))
end end end
result_ieee_cable = HPM.solve_hpf(hdata_ieee_cable, PMs.IVRPowerModel, solver_nlp)














################################################################
###########      PLOTTING    ###################################
################################################################
################################################################


nb_mv = 17
nb_hv = 4


#### Plot harmonic impedance

zh = plot(1:491, abs.(Zh[nb_mv]), yaxis= :log10, ylim = (0.01, 10),  label = "Case A: Inductive network, MV bus", c = :orange,
xlabel="harmonic h [-]",
ylabel = "\$ Z_{h}\$ in pu", fontfamily = "Computer Modern", xticks = ([1,91,191,291,391,491],["1" "10" "20" "30" "40" "50"])) #, xticklabel = ["$i" for i in 1:50])

plot!(zh, 1:491, abs.(Zh_cable[nb_mv]), yaxis= :log10, ylim = (0.01, 10), label = "Case B: Underground cable network, MV bus", c = :blue,
xlabel="harmonic h [-]",
ylabel = "\$ Z_{h}\$ in pu", fontfamily = "Computer Modern", legend=:bottomright, xticks = ([1,91,191,291,391,491],["1" "10" "20" "30" "40" "50"]))

plot!(zh, 1:491, abs.(Zh[nb_hv]), yaxis= :log10, ylim = (0.01, 10), label = "Case A: Inductive network, HV bus", xticks = ([1,91,191,291,391,491],["1" "10" "20" "30" "40" "50"]),  c = :orange, linestyle = :dash,
xlabel="harmonic h [-]",
ylabel = "\$ Z_{h}\$ in pu", fontfamily = "Computer Modern")

plot!(zh, 1:491, abs.(Zh_cable[nb_hv]), yaxis= :log10, ylim = (0.01, 10), label = "Case B: Underground cable network, HV bus", xticks = ([1,91,191,291,391,491],["1" "10" "20" "30" "40" "50"]),  c = :blue, linestyle = :dash,
xlabel="harmonic h [-]",
ylabel = "\$ Z_{h}\$ in pu", fontfamily = "Computer Modern", legend=:bottomright, yminorgrid = true)
savefig(zh, joinpath(HPM.BASE_DIR,"results/zh.pdf"))


### Harmonic voltages
uh_ieee_mv = [vm(result_ieee["solution"]["nw"]["$nh"]["bus"]["$nb_mv"]) for nh in 2:50]
uh_iec_mv = [vm(result_iec["solution"]["nw"]["$nh"]["bus"]["$nb_mv"]) for nh in 2:50]
uh_hhc_mv = [vm(result_hhc["solution"]["nw"]["$nh"]["bus"]["$nb_mv"]) for nh in 2:50]

uh_ieee_mv_c = [vm(result_ieee_cable["solution"]["nw"]["$nh"]["bus"]["$nb_mv"]) for nh in 2:50]
uh_iec_mv_c = [vm(result_iec_cable["solution"]["nw"]["$nh"]["bus"]["$nb_mv"]) for nh in 2:50]
uh_hhc_mv_c = [vm(result_hhc_cable["solution"]["nw"]["$nh"]["bus"]["$nb_mv"]) for nh in 2:50]


uh_ieee_hv = [vm(result_ieee["solution"]["nw"]["$nh"]["bus"]["$nb_hv"]) for nh in 2:50]
uh_iec_hv = [vm(result_iec["solution"]["nw"]["$nh"]["bus"]["$nb_hv"]) for nh in 2:50]
uh_hhc_hv = [vm(result_hhc["solution"]["nw"]["$nh"]["bus"]["$nb_hv"]) for nh in 2:50]

uh_ieee_hv_c = [vm(result_ieee_cable["solution"]["nw"]["$nh"]["bus"]["$nb_hv"]) for nh in 2:50]
uh_iec_hv_c = [vm(result_iec_cable["solution"]["nw"]["$nh"]["bus"]["$nb_hv"]) for nh in 2:50]
uh_hhc_hv_c = [vm(result_hhc_cable["solution"]["nw"]["$nh"]["bus"]["$nb_hv"]) for nh in 2:50]

uh_ind = plot([0.0, 52.0], [0.03 0.03]', label = "Limit MV", linewidth = 3, c = :black)
plot!(uh_ind, [0.0, 52.0], [0.015 0.015]', label = "Limit HV", linewidth = 3, c = :black, linestyle = :dot)
plot!(uh_ind, 2:50, uh_ieee_mv,  marker = :diamond, c = :orange,
xlabel="harmonic h [-]",
ylabel="\$𝗨^{ihd}\$ [pu]",
fontfamily = "Computer Modern",
size = (1200, 266),
left_margin = 10mm, bottom_margin = 10mm,
ylim = (0.0, 0.08),
label = "IEEE519-2022 MV"
)

plot!(uh_ind, 2:50, uh_iec_mv,  marker = :diamond, c = :blue,
xlabel="harmonic h [-]",
ylabel="\$𝗨^{ihd}\$ [pu]",
fontfamily = "Computer Modern",
size = (1200, 266),
left_margin = 10mm, bottom_margin = 10mm,
ylim = (0.0, 0.08),
label = "IEC61000-3-6:2008 MV"
)

plot!(uh_ind, 2:50, uh_hhc_mv,  marker = :diamond, c = :green,
xlabel="harmonic h [-]",
ylabel="\$𝗨^{ihd}\$ [pu]",
fontfamily = "Computer Modern",
size = (1200, 266),
left_margin = 10mm, bottom_margin = 10mm,
ylim = (0.0, 0.08),
label = "Optimisation approach MV"
)

plot!(uh_ind, 2:50, uh_ieee_hv,  marker = :circle, c = :orange, linestyle = :dot,
xlabel="harmonic h [-]",
ylabel="\$𝗨^{ihd}\$ [pu]",
fontfamily = "Computer Modern",
size = (1200, 266),
left_margin = 10mm, bottom_margin = 10mm,
ylim = (0.0, 0.08),
label = "IEEE519-2022 HV"
)


plot!(uh_ind, 2:50, uh_iec_hv,  marker = :circle, c = :blue, linestyle = :dot,
xlabel="harmonic h [-]",
ylabel="\$𝗨^{ihd}\$ [pu]",
fontfamily = "Computer Modern",
size = (1200, 266),
left_margin = 10mm, bottom_margin = 10mm,
ylim = (0.0, 0.08),
label = "IEC61000-3-6:2008 HV"
)


plot!(uh_ind, 2:50, uh_hhc_hv,  marker = :circle, c = :green, linestyle = :dot,
xlabel="harmonic h [-]",
ylabel="\$𝗨^{ihd}\$ [pu]",
fontfamily = "Computer Modern",
size = (1200, 266),
left_margin = 10mm, bottom_margin = 10mm,
ylim = (0.001, 0.1),
label = "Optimisation approach HV", yaxis = :log10, legend = :outertopright, yminorgrid=true
)

savefig(uh_ind, joinpath(HPM.BASE_DIR,"results/uhd_ind.pdf"))

uh_cap = plot([0.0, 52.0], [0.03 0.03]', label = "Limit MV", linewidth = 3, c = :black)
plot!(uh_cap, [0.0, 52.0], [0.015 0.015]', label = "Limit HV", linewidth = 3, c = :black, linestyle = :dot)
plot!(uh_cap,2:50, uh_ieee_mv_c,  marker = :diamond, xlabel="harmonic h [-]", c = :orange,
ylabel="\$𝗨^{ihd}\$ [pu]",
fontfamily = "Computer Modern",
size = (1200, 266),
left_margin = 10mm, bottom_margin = 10mm,
ylim = (0.0, 0.08),
label = "IEEE519-2022 MV"
)

plot!(uh_cap, 2:50, uh_iec_mv_c,  marker = :diamond, xlabel="harmonic h [-]", c = :blue,
ylabel="\$𝗨^{ihd}\$ [pu]",
fontfamily = "Computer Modern",
size = (1200, 266),
left_margin = 10mm, bottom_margin = 10mm,
ylim = (0.0, 0.08),
label = "IEC61000-3-6:2008 MV"
)

plot!(uh_cap, 2:50, uh_hhc_mv_c,  marker = :diamond, xlabel="harmonic h [-]", c = :green,
ylabel="\$𝗨^{ihd}\$ [pu]",
fontfamily = "Computer Modern",
size = (1200, 266),
left_margin = 10mm, bottom_margin = 10mm,
ylim = (0.0, 0.08),
label = "Optimisation approach MV"
)

plot!(uh_cap, 2:50, uh_ieee_hv_c,  marker = :circle, xlabel="harmonic h [-]", c = :orange, linestyle = :dot,
ylabel="\$𝗨^{ihd}\$ [pu]",
fontfamily = "Computer Modern",
size = (1200, 266),
left_margin = 10mm, bottom_margin = 10mm,
ylim = (0.0, 0.08),
label = "IEEE519-2022 HV"
)


plot!(uh_cap, 2:50, uh_iec_hv_c,  marker = :circle, xlabel="harmonic h [-]", c = :blue, linestyle = :dot,
ylabel="\$𝗨^{ihd}\$ [pu]",
fontfamily = "Computer Modern",
size = (1200, 266),
left_margin = 10mm, bottom_margin = 10mm,
ylim = (0.0, 0.08),
label = "IEC61000-3-6:2008 HV"
)


plot!(uh_cap, 2:50, uh_hhc_hv_c,  marker = :circle, xlabel="harmonic h [-]", c = :green, linestyle = :dot,
ylabel="\$𝗨^{ihd}\$ [pu]",
fontfamily = "Computer Modern",
size = (1200, 266),
left_margin = 10mm, bottom_margin = 10mm,
ylim = (0.0001, 0.2),
label = "Optimisation approach HV", yaxis = :log10, legend = :outertopright, yminorgrid=true
)

savefig(uh_cap, joinpath(HPM.BASE_DIR,"results/uhd_cap.pdf"))




###### Harmonic currents
mv_load = 12
hv_load = 3

ih_ieee_mv = [result_ieee["solution"]["nw"]["$nh"]["load"]["$mv_load"]["cm"] for nh in 2:50]
ih_iec_mv = [result_iec["solution"]["nw"]["$nh"]["load"]["$mv_load"]["cm"] for nh in 2:50]
ih_hhc_mv = [result_hhc["solution"]["nw"]["$nh"]["load"]["$mv_load"]["cm"] for nh in 2:50]

ih_ieee_mv_c = [result_ieee_cable["solution"]["nw"]["$nh"]["load"]["$mv_load"]["cm"] for nh in 2:50]
ih_iec_mv_c = [result_iec_cable["solution"]["nw"]["$nh"]["load"]["$mv_load"]["cm"] for nh in 2:50]
ih_hhc_mv_c = [result_hhc_cable["solution"]["nw"]["$nh"]["load"]["$mv_load"]["cm"] for nh in 2:50]

ih_ieee_hv = [result_ieee["solution"]["nw"]["$nh"]["load"]["$hv_load"]["cm"] for nh in 2:50]
ih_iec_hv = [result_iec["solution"]["nw"]["$nh"]["load"]["$hv_load"]["cm"] for nh in 2:50]
ih_hhc_hv = [result_hhc["solution"]["nw"]["$nh"]["load"]["$hv_load"]["cm"] for nh in 2:50]

ih_ieee_hv_c = [result_ieee_cable["solution"]["nw"]["$nh"]["load"]["$hv_load"]["cm"] for nh in 2:50]
ih_iec_hv_c = [result_iec_cable["solution"]["nw"]["$nh"]["load"]["$hv_load"]["cm"] for nh in 2:50]
ih_hhc_hv_c = [result_hhc_cable["solution"]["nw"]["$nh"]["load"]["$hv_load"]["cm"] for nh in 2:50]



ih_ind = plot(2:50, ih_ieee_mv,  marker = :diamond,  c = :orange,
xlabel="harmonic h [-]",
ylabel="\$𝗜^{ihd}\$ [pu]",
fontfamily = "Computer Modern",
size = (1200, 266),
left_margin = 10mm, bottom_margin = 10mm,
ylim = (0.0001, 0.01),
label = "IEEE519-2022 MV"
)

plot!(ih_ind, 2:50, ih_iec_mv,  marker = :diamond,c = :blue,
xlabel="harmonic h [-]", 
ylabel="\$𝗜^{ihd}\$ [pu]",
fontfamily = "Computer Modern",
size = (1200, 266),
left_margin = 10mm, bottom_margin = 10mm,
ylim = (0.0001, 0.01),
label = "IEC61000-3-6:2008 MV"
)

plot!(ih_ind, 2:50, ih_hhc_mv,  marker = :diamond, c = :green,
xlabel="harmonic h [-]",
ylabel="\$𝗜^{ihd}\$ [pu]",
fontfamily = "Computer Modern",
size = (1200, 266),
left_margin = 10mm, bottom_margin = 10mm,
ylim = (0.0001, 0.01),
label = "Optimisation approach MV"
)

plot!(ih_ind, 2:50, ih_ieee_hv,  marker = :circle, c = :orange, linestyle = :dot,
xlabel="harmonic h [-]",
ylabel="\$𝗜^{ihd}\$ [pu]",
fontfamily = "Computer Modern",
size = (1200, 266),
left_margin = 10mm, bottom_margin = 10mm,
ylim = (0.0001, 0.01),
label = "IEEE519-2022 HV"
)


plot!(ih_ind, 2:50, ih_iec_hv,  marker = :circle, c = :blue, linestyle = :dot,
xlabel="harmonic h [-]", 
ylabel="\$𝗜^{ihd}\$ [pu]",
fontfamily = "Computer Modern",
size = (1200, 266),
left_margin = 10mm, bottom_margin = 10mm,
ylim = (0.0001, 0.01),
label = "IEC61000-3-6:2008 HV"
)


plot!(ih_ind, 2:50, ih_hhc_hv,  marker = :circle,  c = :green, linestyle = :dot,
xlabel="harmonic h [-]",
ylabel="\$𝗜^{ihd}\$ [pu]",
fontfamily = "Computer Modern",
size = (1200, 266),
left_margin = 10mm, bottom_margin = 10mm,
ylim = (0.0001, 0.01),
label = "Optimisation approach HV", yaxis = :log10, legend = :outertopright, yminorgrid=true
)

savefig(ih_ind, joinpath(HPM.BASE_DIR,"results/ihd_ind.pdf"))


ih_cap = plot(2:50, ih_ieee_mv_c,  marker = :diamond, c = :orange,
xlabel="harmonic h [-]",
ylabel="\$𝗜^{ihd}\$ [pu]",
fontfamily = "Computer Modern",
size = (1200, 266),
left_margin = 10mm, bottom_margin = 10mm,
ylim = (0.00001, 0.01),
label = "IEEE519-2022 MV"
)

plot!(ih_cap, 2:50, ih_iec_mv_c,  marker = :diamond,  c = :blue,
xlabel="harmonic h [-]",
ylabel="\$𝗜^{ihd}\$ [pu]",
fontfamily = "Computer Modern",
size = (1200, 266),
left_margin = 10mm, bottom_margin = 10mm,
ylim = (0.00001, 0.01),
label = "IEC61000-3-6:2008 MV"
)

plot!(ih_cap, 2:50, ih_hhc_mv_c,  marker = :diamond,  c = :green,
xlabel="harmonic h [-]",
ylabel="\$𝗜^{ihd}\$ [pu]",
fontfamily = "Computer Modern",
size = (1200, 266),
left_margin = 10mm, bottom_margin = 10mm,
ylim = (0.00001, 0.01),
label = "Optimisation approach MV"
)

plot!(ih_cap, 2:50, ih_ieee_hv_c,  marker = :circle,c = :orange, linestyle = :dot,
xlabel="harmonic h [-]", 
ylabel="\$𝗜^{ihd}\$ [pu]",
fontfamily = "Computer Modern",
size = (1200, 266),
left_margin = 10mm, bottom_margin = 10mm,
ylim = (0.00001, 0.01),
label = "IEEE519-2022 HV"
)


plot!(ih_cap, 2:50, ih_iec_hv_c,  marker = :circle, c = :blue, linestyle = :dot,
xlabel="harmonic h [-]",
ylabel="\$𝗜^{ihd}\$ [pu]",
fontfamily = "Computer Modern",
size = (1200, 266),
left_margin = 10mm, bottom_margin = 10mm,
ylim = (0.00001, 0.01),
label = "IEC61000-3-6:2008 HV"
)


plot!(ih_cap, 2:50, ih_hhc_hv_c,  marker = :circle, c = :green, linestyle = :dot,
xlabel="harmonic h [-]", 
ylabel="\$𝗜^{ihd}\$ [pu]",
fontfamily = "Computer Modern",
size = (1200, 266),
left_margin = 10mm, bottom_margin = 10mm,
ylim = (0.00001, 0.01),
label = "Optimisation approach HV", yaxis = :log10, legend = :outertopright, yminorgrid=true
)

savefig(ih_cap, joinpath(HPM.BASE_DIR,"results/ihd_cap.pdf"))





















# uh = [uh_ieee uh_iec uh_hhc_c uh_ieee_c uh_iec uh_hhc_c]
# ctg = repeat(["IEEE519-2022 Ind." ,"IEC61000-3-6:2008 Ind.", "Optimisation approach Ind.", "IEEE519-2022 Orig." ,"IEC61000-3-6:2008 Orig.", "Optimisation approach Orig."], inner = 49)
# nam = repeat(collect(2:50), outer = 6)
# p_hv = groupedbar(
#     nam, uh, bar_position = :dodge, bar_width = 0.7, group = ctg,
#     lw = 0, framestyle = :box,
#     xlabel="harmonic h [-]",
#     ylabel="\$𝗨^{ihd}\$ [pu]",
#     fontfamily = "Computer Modern",
#     size = (1200, 266),
#     left_margin = 10mm, bottom_margin = 10mm,
#     ylim = (0.0, 0.15)
# )

# savefig(p_hv, joinpath(HPM.BASE_DIR,"results/uhd_hv.pdf"))


# uh = [uh_ieee uh_iec uh_hhc_c uh_ieee_c uh_iec uh_hhc_c]
# ctg = repeat(["IEEE519-2022 Ind." ,"IEC61000-3-6:2008 Ind.", "Optimisation approach Ind.", "IEEE519-2022 Orig." ,"IEC61000-3-6:2008 Orig.", "Optimisation approach Orig."], inner = 49)
# nam = repeat(collect(2:50), outer = 6)
# p_mw = groupedbar(
#     nam, uh, bar_position = :dodge, bar_width = 0.7, group = ctg,
#     lw = 0, framestyle = :box,
#     xlabel="harmonic h [-]",
#     ylabel="\$𝗨^{ihd}\$ [pu]",
#     fontfamily = "Computer Modern",
#     size = (1200, 266),
#     left_margin = 10mm, bottom_margin = 10mm,
#     ylim = (0.0, 0.08)
# )

# savefig(p_mw, joinpath(HPM.BASE_DIR,"results/uhd_mv.pdf"))

# ih = [ih_ieee ih_iec ih_hhc ih_ieee_c ih_iec_c ih_hhc_c]
# ctg = repeat(["IEEE519-2022 Ind." ,"IEC61000-3-6:2008 Ind.", "Optimisation approach Ind.", "IEEE519-2022 Orig." ,"IEC61000-3-6:2008 Orig.", "Optimisation approach Orig."], inner = 49)
# nam = repeat(collect(2:50), outer = 6)
# pi_mv = groupedbar(
#     nam, ih, bar_position = :dodge, bar_width = 0.7, group = ctg,
#     lw = 0, framestyle = :box,
#     xlabel="harmonic h [-]",
#     ylabel="\$𝗜^{ihd}\$ [pu]",
#     fontfamily = "Computer Modern",
#     size = (1200, 266),
#     left_margin = 10mm, bottom_margin = 10mm,
#     ylim = (0.0, 0.01)
# )

# savefig(pi_mv, joinpath(HPM.BASE_DIR,"results/ihd_mv.pdf"))




# ih = [ih_ieee ih_iec ih_hhc ih_ieee_c ih_iec_c ih_hhc_c]
# ctg = repeat(["IEEE519-2022 Ind." ,"IEC61000-3-6:2008 Ind.", "Optimisation approach Ind.", "IEEE519-2022 Orig." ,"IEC61000-3-6:2008 Orig.", "Optimisation approach Orig."], inner = 49)
# nam = repeat(collect(2:50), outer = 6)
# pi_hv = groupedbar(
#     nam, ih, bar_position = :dodge, bar_width = 0.7, group = ctg,
#     lw = 0, framestyle = :box,
#     xlabel="harmonic h [-]",
#     ylabel="\$𝗜^{ihd}\$ [pu]",
#     fontfamily = "Computer Modern",
#     size = (1200, 266),
#     left_margin = 10mm, bottom_margin = 10mm,
#     ylim = (0.0, 0.01)
# )

# savefig(pi_hv, joinpath(HPM.BASE_DIR,"results/ihd_hv.pdf"))










###### Calculate harmonic currents based on the standard:
# -> each harmonic voltage needs to be less than 3% of fundamental
# -> THD below 5%


# include(joinpath(HPM.BASE_DIR,"src/util/hi.jl"))
# path = joinpath(HPM.BASE_DIR,"test/data/matpower/30bus_network_hhc.m")

# path = joinpath(HPM.BASE_DIR,"test/data/matpower/30bus_network_hhc_519.m")
# data = PMs.parse_file(path)

# data_cable = deepcopy(data)

# ihd_limits = [1.00000, 0.02000, 0.05000, 0.01000, 0.06000, 
# 0.00500, 0.05000, 0.00500, 0.01500, 0.00500, 
# 0.03500, 0.00458, 0.03000, 0.00429, 0.00400, 
# 0.00406, 0.02000, 0.00389, 0.01761, 0.00375, 
# 0.00200, 0.00364, 0.01408, 0.00354, 0.01275, 
# 0.00346, 0.00200, 0.00339, 0.01061, 0.00333, 
# 0.00975, 0.00328, 0.00200, 0.00324, 0.00833, 
# 0.00319, 0.00773, 0.00316, 0.00200, 0.00313, 
# 0.00671, 0.00310, 0.00627, 0.00307, 0.00200, 
# 0.00304, 0.00551, 0.00302, 0.00518, 0.00300]

# # reset all susceptance of the branches
# for (nb, branch) in data["branch"] 
#     branch["b_fr"] = 0.0
#     branch["b_to"] = 0.0
# end

# # add the load current magnitude [pu] to all buses
# for (nb, bus) in data_cable["bus"]
#     bus["i_load"] = 0.0
# end
# for (nl, load) in data_cable["load"]
#     nb = load["load_bus"]
#     data["bus"]["$nb"]["i_load"] = abs(load["pd"] + im * load["qd"])
# end

# # add the short-circuit current magnitude [pu] to all busses
# i_sc = [0.163, 14.94, 11.828, 12.944, 4.849, 10.76, 8.024, 5.19, 6.093, 6.09, 
#         2.366, 5.342, 2.556, 3.55, 4.64, 4.234, 4.713, 3.426, 3.42, 3.669, 4.82, 
#         4.806, 3.516, 4.033, 2.89, 1.312, 3.117, 8.267, 1.601, 1.461]
# for (nb, bus) in data_cable["bus"]
#     bus["i_sc"] = i_sc[parse(Int, nb)]
# end

# # add the short-circuit vs load current ratio
# for (nb, bus) in data_cable["bus"]
#     bus["i_ratio"] = bus["i_sc"] / bus["i_load"]
#     bus["i_ihd_iec"] = [ihd_limits[h] * bus["i_load"] for h in 2:50]
# end

# hdata_ind  = HPM.replicate(data_cable, H=collect(1:50))
# for (nw, ntw) in hdata_ind["nw"] if parse(Int, nw) ≠ 1
#     for (nl, load) in ntw["load"]
#         nb = load["load_bus"]

#         #v = data["bus"]["$nb"]["base_kv"]
#         #i = data["bus"]["$nb"]["i_ratio"]

#         load["multiplier"] = ihd_limits[parse(Int,nw)]
# end end end
# results_hpf_ind = HPM.solve_hpf(hdata_ind, PMs.IVRPowerModel, solver_nlp)




# # Plots
# ## input
# nb = 17
# nl = 1
# ## function


# ## plots 
# # bar(2:50, [data_cable["bus"]["$nb"]["i_ihd_iec"][nh-1] for nh in 2:50], 
# #         label="IEC61000-3-6:2008",
# #         xlabel="harmonic h [-]",
# #         ylabel="individual harmonic unit current Iⁱʰᵈ [pu]")
# # bar(2:50, [results_hpf_ind["solution"]["nw"]["$nh"]["load"]["$nl"]["cm"] for nh in 2:50])
# bar(1:50, [vm(results_hpf_ind["solution"]["nw"]["$nh"]["bus"]["$nb"]) for nh in 2:50],
#         label="IEC61000-3-6:2008",
#         xlabel="harmonic h [-]",
#         ylabel="individual harmonic bus voltage Uⁱʰᵈ [pu]")


# r_hhc =  result_hhc["0"]
# ## plots 
# # bar(2:50, [data_cable["bus"]["$nb"]["i_ihd_iec"][nh-1] for nh in 2:50], 
# #         label="Optimisation approach",
# #         xlabel="harmonic h [-]",
# #         ylabel="individual harmonic unit current Iⁱʰᵈ [pu]")
# # bar(2:50, [results_hpf_ind["solution"]["nw"]["$nh"]["load"]["$nl"]["cm"] for nh in 2:50])
# bar(1:50, [vm(r_hhc["solution"]["nw"]["$nh"]["bus"]["$nb"]) for nh in 2:50],
#         label="Optimisation approach",
#         xlabel="harmonic h [-]",
#         ylabel="individual harmonic bus voltage Uⁱʰᵈ [pu]")


# bh_iec = zeros(50, 30)
# hb_ids = zeros(50,30)
# for nb in sort(collect(parse.(Int, keys(results_hpf_ind["solution"]["nw"]["1"]["bus"]))))
#     hb_ids[:, nb] = 1:50
# end
# for nh in 2:50
#     for nb in sort(collect(parse.(Int, keys(results_hpf_ind["solution"]["nw"]["$nh"]["bus"]))))
#     bh_iec[nh, nb] = vm(results_hpf_ind["solution"]["nw"]["$nh"]["bus"]["$nb"])
#     end
# end

# scatter(1:50,)




# # H = 50:50.0:2500
# # ih = Dict{String, Any}(["$i" => zeros(length(H)) for i in cable_branches])
# # ih_min = Dict{String, Any}(["$i" => zeros(length(H)) for i in cable_branches])
# # ih_lin = Dict{String, Any}(["$i" => zeros(length(H)) for i in cable_branches])
# # ih_min_lin = Dict{String, Any}(["$i" => zeros(length(H)) for i in cable_branches])

# # load_buses = [data["load"]["$l"]["load_bus"] for l in sort(parse.(Int, keys(data["load"])))]
# # global_contribution = [abs.(data["load"]["$l"]["pd"] + im * data["load"]["$l"]["qd"]) for l in sort(parse.(Int, keys(data["load"])))] ./ sum([abs.(load["pd"] + im * load["qd"]) for (l, load) in data["load"]])




# # planning_levels_hv = [
# # 0.01400, 0.00800, 0.02000, 0.00400, 0.02000,
# # 0.00400, 0.02000, 0.00400, 0.01000, 0.00350,
# # 0.01500, 0.00318, 0.01500, 0.00296, 0.00300,
# # 0.00279, 0.01200, 0.00266, 0.01074, 0.00255,
# # 0.00200, 0.00246, 0.00887, 0.00239, 0.00816,
# # 0.00233, 0.00200, 0.00228, 0.00703, 0.00223,
# # 0.00658, 0.00219, 0.00200, 0.00216, 0.00583,
# # 0.00213, 0.00551, 0.00537, 0.00200, 0.00510,
# # 0.00498, 0.00205, 0.00474, 0.00203, 0.00200,
# # 0.00201, 0.00434, 0.00199, 0.00416, 0.00198
# # ]


# # General summation law:
# # I = ∑i  J[i,h]^α)^(1/α):
# # h < 5 -> α = 1 
# # 5 <= h < 10 -> α = 1.4
# # 5 <= h < 10 -> α = 2

# # # Calculation with actual harmonic impedance
# # for c in cable_branches
# #     Zh = Zc["$c"]
# #     for h in 2:length(H)
# #         if h < 5
# #             α = 1
# #         elseif 5 <= h && h <= 10
# #             α = 1.4
# #         else
# #             α = 2
# #         end
# #         α = 1
# #         ih["$c"][h] = sum([abs(ihd_limits[h] / Zh[load_buses[i]][h])^α  * global_contribution[i] for i in 1:length(load_buses)])^(1/α)
# #         ih_min["$c"][h] = minimum([abs(ihd_limits[h] / Zh[load_buses[i]][h]) * global_contribution[i]  for i in 1:length(load_buses)]')
# #     end
# # end

# # # Calculation with linear impdedance assumption harmonic impedance
# # for c in cable_branches
# #     Zh = Zc["$c"]
# #     for h in 2:length(H)
# #         if h < 5
# #             α = 1
# #         elseif 5 <= h && h <= 10
# #             α = 1.4
# #         else
# #             α = 2
# #         end
# #         α = 1
# #         ih_lin["$c"][h] = sum([abs(ihd_limits[h] / (h * Zh[load_buses[i]][1]))^α  * global_contribution[i] for i in 1:length(load_buses)])^(1/α)
# #         ih_min_lin["$c"][h] = minimum([abs(ihd_limits[h] / (h * Zh[load_buses[i]][1])) * global_contribution[i] for i in 1:length(load_buses)]')
# #     end
# # end

# # thd_ih = [sqrt(sum(ih["$c"][2:end].^2)) for c in cable_branches]
# # sum_ih_min = [sum(ih_min["$c"][2:end]) for c in cable_branches] .* length(load_buses)
# # thd_ih_lin = [sqrt(sum(ih_lin["$c"][2:end].^2)) for c in cable_branches]
# # sum_ih_lin = [sum(ih_lin["$c"][2:end]) for c in cable_branches]
# # sum_ih_min_lin = [sum(ih_min_lin["$c"][2:end]) for c in cable_branches] .* length(load_buses)


# # plot!(hhc, Int.(0:length(cable_branches)-1), thd_ih', xlabel = "Number of HV cables", ylabel = "HHC in pu", fontfamily = "Computer Modern", label = "IEC61000-3-6:2008", ylim = (0,1))


# # savefig(hhc, joinpath(HPM.BASE_DIR,"results/hosting_capacity.pdf"))



