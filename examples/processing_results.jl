using Plots
using JSON
using HarmonicPowerModels

const HPM = HarmonicPowerModels

case = "case1888_rte"
harmonic_range = "__1_50_"

# Laoding input data

path = joinpath(HPM.BASE_DIR,"test","data","matpower", join([case, ".m"]))
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

filename = joinpath(HPM.BASE_DIR, "results", join([case, harmonic_range, "Kalai-Smorodinsky bargaining.json"]))
r_ksb = Dict{String, Any}()
open(filename) do f
dicttxt = read(f, String)
global r_ksb = JSON.parse(dicttxt)
end


bus_id = 1

vm_ae = zeros(1, 49)
vm_mm = zeros(1, 49)
vm_me = zeros(1, 49)
vm_ksb = zeros(1, 49)

for idx in 1:49
    h = idx + 1
    vm_ae[idx] = r_ae["solution"]["nw"]["$h"]["bus"]["$bus_id"]["vm"]
    vm_mm[idx] = r_mm["solution"]["nw"]["$h"]["bus"]["$bus_id"]["vm"]
    vm_me[idx] = r_me["solution"]["nw"]["$h"]["bus"]["$bus_id"]["vm"]
    vm_ksb[idx] = r_ksb["solution"]["nw"]["$h"]["bus"]["$bus_id"]["vm"]
end

vmh = plot(1:49, vm_ae', marker = :diamond, label = "absolute equality", yaxis = :log10)
plot!(vmh, 1:49, vm_me', marker = :diamond, label = "maximum efficieny", yaxis = :log10)
plot!(vmh, 1:49, vm_ksb', marker = :diamond, label = "Kalai-Smorodinsky bargaining", yaxis = :log10)
plot!(vmh, 1:49, vm_mm', marker = :diamond, label = "maximin", yaxis = :log10, ylabel = "\$V_{h}~[pu]\$", xlabel = "h [-]", fontfamily = "Computer Modern")












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
