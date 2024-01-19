################################################################################
#  Copyright 2023, Hakan Ergun                                                 #
################################################################################
# HarmonicPowerModels.jl                                                       #
# An extention package of PowerModels(Distribution).jl for Harmonics           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################

# using pkgs
using HarmonicPowerModels, Ipopt, JuMP, PowerModels, Plots, StatsPlots, ElectricalEngineering

output_path = "/Users/hergun/Library/CloudStorage/OneDrive-KULeuven/Projects/HARMONIC/WP5/plots"

# pkgs cte
const PMs = PowerModels
const HPM = HarmonicPowerModels
const EE = ElectricalEngineering

# set the solver
solver = Ipopt.Optimizer

# set the formulation
form = dHHC_NLP

# read-in data 
path = joinpath(HPM.BASE_DIR,"test/data/matpower/industrial_network_hhc.m")
data = PMs.parse_file(path)

# build harmonic data
hdata = HPM.replicate(data, H=[1, 3, 5, 7, 9, 11, 13])

# vector group shift
vector_shift = Dict(1 => 0, 2 => 11/6 * pi, 3 => 11/6 * pi, 4 => 11/6 * pi, 5 => 11/6 * pi, 6 => 5/3 * pi, 7 => 5/3 * pi, 8 => 5/3 * pi)

for H=[1, 3, 5, 7, 9, 11, 13]
    for (l, load) in hdata["nw"]["$H"]["load"]
        load["c_rating"] = 1.0
        bus_id = load["load_bus"]
        load["reference_harmonic_angle"] = 3*pi/4 + vector_shift[bus_id]# rad
        error = (pi/20) * (1 + 1 / 10)
        load["harmonic_angle_range"] = error # rad, symmetric around reference
    end
end
results_hhc = HPM.solve_hhc(hdata, form, solver)
HPM.plot_harmonic_current_injections(results_hhc, 3, output_path)
HPM.plot_harmonic_voltages(results_hhc, 11, output_path)

ang_pos = [0 pi/4 pi/2 3*pi/4 pi 5*pi/4 6*pi/4 7*pi/4]
objective = zeros(10,8)
idx = 1
for p in ang_pos
    for idx_1 in 1:10
        # Define ref_angle and angle range 
        for H=[1, 3, 5, 7, 9, 13]
            for (l, load) in hdata["nw"]["$H"]["load"]
                load["c_rating"] = 1.0
                bus_id = load["load_bus"]
                load["reference_harmonic_angle"] = p + vector_shift[bus_id]# rad
                error = (pi/20) * (1 + idx_1 / 10)
                load["harmonic_angle_range"] = error # rad, symmetric around reference
            end
        end
        # solve HC problem
        results_hhc = HPM.solve_hhc(hdata, form, solver)
        objective[idx_1, idx] = results_hhc["objective"]
    end
    global idx = idx + 1
end

a = zeros(length(objective),1)
a = objective[:,1]

p1 = Plots.scatter(ang_pos, objective, legend = nothing, xlabel = "\$Harmonic~injection~angle\$", ylabel = "\$Total~HHC\$", xticks = (ang_pos, ["0", "π/4", "π/2", "3π/4", "π", "5π/4", "3π/2", "7π/4"]), xtickfont = "Computer Modern", ytickfont = "Computer Modern", fontfamily = "Computer Modern")
Plots.savefig(p1, joinpath(output_path,"angle_vs_hhc.pdf"))

p2 = StatsPlots.boxplot(ang_pos, objective, legend = nothing, xlabel = "\$Harmonic~injection~angle\$", ylabel = "\$Total~HHC\$", xticks = (ang_pos, ["0", "π/4", "π/2", "3π/4", "π", "5π/4", "3π/2", "7π/4"]), xtickfont = "Computer Modern", ytickfont = "Computer Modern", fontfamily = "Computer Modern")
Plots.savefig(p2, joinpath(output_path,"angle_range_vs_hhc.pdf"))