using Plots
using JSON
using HarmonicPowerModels

const HPM = HarmonicPowerModels

case = "nem2300harmonic"
harmonic_range = "__1_50_"

# Laoding input data

path = joinpath(HPM.BASE_DIR,"test","data","json", join([case, ".json"]))
data = Dict{String, Any}()
open(path) do f
dicttxt = read(f,String)  # file information to string
    global data = JSON.parse(dicttxt)  # parse and transform data
end

##### Processing output #######
filename = joinpath(HPM.BASE_DIR, "results", join([case, harmonic_range, "maximum efficiency.json"]))
r_me = Dict{String, Any}()
open(filename) do f
dicttxt = read(f, String)
global r_me = JSON.parse(dicttxt)
end

filename = joinpath(HPM.BASE_DIR, "results", join([case,harmonic_range, "absolute equality.json"]))
r_ae = Dict{String, Any}()
open(filename) do f
dicttxt = read(f, String)
global r_ae = JSON.parse(dicttxt)
end

filename = joinpath(HPM.BASE_DIR, "results", join([case, harmonic_range, "maximin.json"]))
r_mm = Dict{String, Any}()
open(filename) do f
dicttxt = read(f, String)
global r_mm = JSON.parse(dicttxt)
end


bus_id = 













# bus_voltage_magnitudes = zeros(length(H), length(data["bus"]))
# bus_voltage_angles = zeros(length(H), length(data["bus"]))
# harmonic_current_injections = zeros(length(H), length(data["load"]))

# for h_idx in 2:length(H)
#     h = H[h_idx]
#     for (b, bus) in r_ae["solution"]["nw"]["$h"]["bus"]
#         bus_voltage_magnitudes[h_idx, parse(Int, b)] = bus["vm"]
#         bus_voltage_angles[h_idx, parse(Int, b)] = bus["va"]
#     end
#     for (l, load) in r_ae["solution"]["nw"]["$h"]["load"]
#         harmonic_current_injections[h_idx, parse(Int, l)] = load["cm"]
#     end
# end

# StatsPlots.boxplot(bus_voltage_magnitudes')
# StatsPlots.boxplot(bus_voltage_angles')
# StatsPlots.boxplot(harmonic_current_injections', outliers = false)
