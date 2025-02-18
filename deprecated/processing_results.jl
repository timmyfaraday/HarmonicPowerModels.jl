using Plots
using JSON
using HarmonicPowerModels
using PowerModels
using Statistics
using DelimitedFiles


const PMs = PowerModels
const HPM = HarmonicPowerModels


case = "case1888_rte"
harmonic_range = "__1_50_"

# Laoding input data

path = joinpath(HPM.BASE_DIR,"test","data","matpower", join([case, ".m"]))
data = PMs.parse_file(path)
for (b, bus) in data["bus"]
    bus["standard"] = "IEC61000-3-6:2008"
    bus["ref_angle"] = 0.0
end

idx = 1
for b in collect(sort(parse.(Int, keys(data["bus"]))))
    data["bus"]["$b"]["hb_idx"] = idx
    global idx = idx + 1
end

# Enforce bus voltage limits and standard
for (b, bus) in data["bus"]
    bus["standard"] = "IEC61000-3-6:2008"
    bus["ref_angle"] = 0.0
    bus["vmax"] = 1.1#min(1.1, bus["vmax"])
    bus["vmin"] = 0.9##max(0.90, bus["vmin"]) 
end

# get rid of negative reactances and resistances
for (br, branch) in data["branch"]
    branch["tap"] = 1.0
    branch["shift"] = 0.0
    branch["c_rating"] = (branch["rate_a"] / 0.9) 
end

# Switch on all generators and set their lower bound to zero
for (g, gen) in data["gen"]
    gen[ "gen_status"] = 1
    gen["pmin"] = 0.0
    gen_bus = gen["gen_bus"]
    data["bus"]["$gen_bus"]["bus_type"] = 2
end

idx = 1
for b in collect(sort(parse.(Int, keys(data["bus"]))))
    data["bus"]["$b"]["hb_idx"] = idx
    global idx = idx + 1
end

# define the set of considered harmonics
H = [h for h in 1:50]

# solve HHC problem -- SOC 
hdata = HPM.replicate(data, H=H)



##### Processing output #######
filename = joinpath(HPM.BASE_DIR, "results", join([case, "_ihdonly", harmonic_range, "maximum efficiency.json"]))
r_me = Dict{String, Any}()
open(filename) do f
dicttxt = read(f, String)
global r_me = JSON.parse(dicttxt)
end

filename = joinpath(HPM.BASE_DIR, "results", join([case, "_ihdonly",harmonic_range, "absolute equality.json"]))
r_ae = Dict{String, Any}()
open(filename) do f
dicttxt = read(f, String)
global r_ae = JSON.parse(dicttxt)
end

filename = joinpath(HPM.BASE_DIR, "results", join([case, "_ihdonly", harmonic_range, "maximin.json"]))
r_mm = Dict{String, Any}()
open(filename) do f
dicttxt = read(f, String)
global r_mm = JSON.parse(dicttxt)
end

filename = joinpath(HPM.BASE_DIR, "results", join([case, "_ihdonly", harmonic_range, "Kalai-Smorodinsky bargaining.json"]))
r_ksb = Dict{String, Any}()
open(filename) do f
dicttxt = read(f, String)
global r_ksb = JSON.parse(dicttxt)
end


##### Harmonic current injections per node:


load_gen_buses = []

filter_gens = false

ih_me = zeros(50, length(r_me["solution"]["nw"]["2"]["load"]))
ih_ae = zeros(50, length(r_me["solution"]["nw"]["2"]["load"]))
ih_mm = zeros(50, length(r_me["solution"]["nw"]["2"]["load"]))
ih_ksb = zeros(50, length(r_ksb["solution"]["nw"]["2"]["load"]))


for (l, load) in data["load"]
    for (g, gen) in data["gen"]
        if gen["gen_bus"] == load["load_bus"]
            push!(load_gen_buses,  gen["gen_bus"])
        end
    end
end

for (n, nw) in r_me["solution"]["nw"]
    if n ≠ "1"
        for (l, load) in nw["load"]
            if !any(data["load"][l]["load_bus"] .== load_gen_buses) || filter_gens == false
                l_ = parse(Int, l)
                n_ = parse(Int, n)
                ih_me[n_, l_] = load["cm"]
            end
        end
    end
end

for (n, nw) in r_ae["solution"]["nw"]
    if n ≠ "1"
        for (l, load) in nw["load"]
            if !any(data["load"][l]["load_bus"] .== load_gen_buses) || filter_gens == false
                l_ = parse(Int, l)
                n_ = parse(Int, n)
                ih_ae[n_, l_] = load["cm"]
            end
        end
    end
end

for (n, nw) in r_mm["solution"]["nw"]
    if n ≠ "1"
        for (l, load) in nw["load"]
            if !any(data["load"][l]["load_bus"] .== load_gen_buses) || filter_gens == false
                l_ = parse(Int, l)
                n_ = parse(Int, n)
                ih_mm[n_, l_] = load["cm"]
            end
        end
    end
end

for (n, nw) in r_ksb["solution"]["nw"]
    if n ≠ "1"
        for (l, load) in nw["load"]
            if !any(data["load"][l]["load_bus"] .== load_gen_buses) || filter_gens == false
                l_ = parse(Int, l)
                n_ = parse(Int, n)
                ih_ksb[n_, l_] = load["cm"]
            end
        end
    end
end

hhc_me = sum(sum(ih_me))
hhc_ae = sum(sum(ih_ae))
hhc_mm = sum(sum(ih_mm))
hhc_ksb = sum(sum(ih_ksb))

println("HHC maximum efficiency: ", hhc_me)
println("HHC maximin: ", hhc_mm)
println("HHC absolute equality: ", hhc_ae)
println("HHC bargaining: ", hhc_ksb)

pih_me = heatmap(log10.(ih_me)', c=cgrad(:jet, rev = true), zlabel = "\$I_{h} in pu \$", ylabel = "Load ID", xlabel = "h [-]", fontfamily = "Computer Modern")
pih_mm = heatmap(log10.(ih_mm)', c=cgrad(:jet, rev = true), zlabel = "\$I_{h} in pu \$", ylabel = "Load ID", xlabel = "h [-]", fontfamily = "Computer Modern")
pih_ae = heatmap(log10.(ih_ae)', c=cgrad(:jet, rev = true), zlabel = "\$I_{h} in pu \$", ylabel = "Load ID", xlabel = "h [-]", fontfamily = "Computer Modern")
pih_ksb = heatmap(log10.(ih_ksb)', c=cgrad(:jet, rev = true), zlabel = "\$I_{h} in pu \$", ylabel = "Load ID", xlabel = "h [-]", fontfamily = "Computer Modern")




# zticks = ([-1 -2 -3 -4 -5 -6], ["\$10^{-1}\$","\$10^{-1}\$","\$10^{-1}\$","\$10^{-1}\$","\$10^{-1}\$","\$10^{-1}\$"])
a = rand(1000, 50)
plot(1:50,1:1000, a, linetype = :heatmap, c=cgrad(:jet, rev = true))

outfile_ksb = joinpath(HPM.BASE_DIR,"results", join(["ih_ksb_no_rms.txt"]))
outfile_ae = joinpath(HPM.BASE_DIR,"results", join(["ih_ae_no_rms.txt"]))
outfile_mm = joinpath(HPM.BASE_DIR,"results", join(["ih_mm_no_rms.txt"]))
outfile_me = joinpath(HPM.BASE_DIR,"results", join(["ih_me_no_rms.txt"]))

f_ksb = open(outfile_ksb, "w")
f_ae = open(outfile_ae, "w")
f_mm = open(outfile_mm, "w")
f_me = open(outfile_me, "w")
for h in 2 : 50
    for u in 1 : 40 : 1000
        println(f_ksb, h, " ", u, " ", ih_ksb[h,u])
        println(f_ae, h, " ", u, " ", ih_ae[h,u])
        println(f_mm, h, " ", u, " ", ih_mm[h,u])
        println(f_me, h, " ", u, " ", ih_me[h,u])
    end
    println(f_ksb)
    println(f_ae)
    println(f_mm)
    println(f_me)
end
close(f_ksb)
close(f_ae)
close(f_mm)
close(f_me)

if filter_gens == true
    savefig(pih_me, joinpath(HPM.BASE_DIR,"results", join(["ihonly_heat_me_filtered.pdf"])))
    savefig(pih_mm, joinpath(HPM.BASE_DIR,"results", join(["ihonly_heat_mm_filtered.pdf"])))
    savefig(pih_ae, joinpath(HPM.BASE_DIR,"results", join(["ihonly_heat_ae_filtered.pdf"])))
    savefig(pih_ksb, joinpath(HPM.BASE_DIR,"results", join(["ihonly_heat_ksb_filtered.pdf"])))
else
    savefig(pih_me, joinpath(HPM.BASE_DIR,"results", join(["ihonly_heat_me.pdf"])))
    savefig(pih_mm, joinpath(HPM.BASE_DIR,"results", join(["ihonly_heat_mm.pdf"])))
    savefig(pih_ae, joinpath(HPM.BASE_DIR,"results", join(["ihonly_heat_ae.pdf"])))
    savefig(pih_ksb, joinpath(HPM.BASE_DIR,"results", join(["ihonly_heat_ksb.pdf"])))
end

######
vm_me = zeros(50, length(r_me["solution"]["nw"]["2"]["bus"]))
vm_ae = zeros(50, length(r_me["solution"]["nw"]["2"]["bus"]))
vm_mm = zeros(50, length(r_me["solution"]["nw"]["2"]["bus"]))
vm_ksb = zeros(50, length(r_me["solution"]["nw"]["2"]["bus"]))
vm_lim = zeros(50, length(r_me["solution"]["nw"]["2"]["bus"]))


for (n, nw) in r_me["solution"]["nw"]
    if n ≠ "1"
        for (b, bus) in nw["bus"]
            b_ = data["bus"][b]["hb_idx"]
            n_ = parse(Int, n)
            vm_me[n_, b_] = bus["vm"]
        end
    end
end

for (n, nw) in r_mm["solution"]["nw"]
    if n ≠ "1"
        for (b, bus) in nw["bus"]
            b_ = data["bus"][b]["hb_idx"]
            n_ = parse(Int, n)
            vm_mm[n_, b_] = bus["vm"]
        end
    end
end

for (n, nw) in r_ae["solution"]["nw"]
    if n ≠ "1"
        for (b, bus) in nw["bus"]
            b_ = data["bus"][b]["hb_idx"]
            n_ = parse(Int, n)
            vm_ae[n_, b_] = bus["vm"]
        end
    end
end

for (n, nw) in r_ksb["solution"]["nw"]
    if n ≠ "1"
        for (b, bus) in nw["bus"]
            b_ = data["bus"][b]["hb_idx"]
            n_ = parse(Int, n)
            vm_ksb[n_, b_] = bus["vm"]
        end
    end
end

for (n, network) in hdata["nw"] 
    if n ≠ "1"
        for (b, bus) in network["bus"]
            b_ = data["bus"][b]["hb_idx"]
            vm_lim[parse(Int, n), b_] = bus["ihdmax"]
        end
    end
end

index = ones(49, 1)
for i in 2:50
    index[i-1] = i
end

avg_vm_ae = Statistics.mean(vm_ae, dims = 2)
avg_vm_mm = Statistics.mean(vm_mm, dims = 2)
avg_vm_me = Statistics.mean(vm_me, dims = 2)
avg_vm_ksb = Statistics.mean(vm_ksb, dims = 2)
avg_vm_lim = Statistics.mean(vm_lim, dims = 2)


avg_vmh = plot(1:49, avg_vm_ae[2:end], marker = :diamond, label = "absolute equality", yaxis = :log10)
plot!(avg_vmh, 1:49, avg_vm_me[2:end], marker = :diamond, label = "maximum efficieny", yaxis = :log10)
plot!(avg_vmh, 1:49, avg_vm_ksb[2:end], marker = :diamond, label = "KS bargaining", yaxis = :log10)
plot!(avg_vmh, 1:49, avg_vm_mm[2:end], marker = :diamond, label = "maximin", yaxis = :log10)
plot!(avg_vmh, 1:49, avg_vm_lim[2:end], linestyle=:dot, marker = :circle, markersize = 2, color = :black, label = "Avg. IEC61000-3-6:2008 limit", yaxis = :log10, 
ylabel = "\$V_{h}~[pu]\$", xlabel = "h [-]", fontfamily = "Computer Modern", legend=:bottomleft)
savefig(avg_vmh, joinpath(HPM.BASE_DIR,"results", join(["avg_vmh.pdf"])))



avg_vm = hcat(index, avg_vm_ae[2:50], avg_vm_me[2:50], avg_vm_mm[2:50], avg_vm_ksb[2:50], avg_vm_lim[2:50])

header = ["h" "ae" "me" "mm" "ksb" "lim"]
writedlm(joinpath(HPM.BASE_DIR,"results", join(["avg_vmh.csv"])),  [header ; avg_vm], ',')



######

big_buses = []


for idx in 1:size(ih_me, 2)
    if ih_me[5, idx] > 0.85
        load_bus = data["load"]["$idx"]["load_bus"]
        voltage = data["bus"]["$load_bus"]["base_kv"]
        gen_bus = false
        for (g, gen) in data["gen"]
            if gen["gen_bus"] == load_bus
                gen_bus = true
            end
        end
        println("load: ", idx , " connected at bus: ", load_bus, " voltage: ", voltage, " generator: ", gen_bus)
        push!(big_buses, load_bus)
    end
end


bus_id = 127

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

vhlim = zeros(1, 50)
for (nw, network) in hdata["nw"] 
    vhlim[parse(Int, nw)] = network["bus"]["$bus_id"]["ihdmax"]
end

vmh = plot(1:49, vm_ae', marker = :diamond, label = "absolute equality", yaxis = :log10)
plot!(vmh, 1:49, vm_me', marker = :diamond, label = "maximum efficieny", yaxis = :log10)
plot!(vmh, 1:49, vm_ksb', marker = :diamond, label = "KS bargaining", yaxis = :log10)
plot!(vmh, 1:49, vm_mm', marker = :diamond, label = "maximin", yaxis = :log10)
plot!(vmh, 1:49, vhlim[2:end], linestyle=:dot, marker = :circle, markersize = 2, color = :black, label = "IEC61000-3-6:2008 limit", yaxis = :log10, 
ylabel = "\$V_{h}~[pu]\$", xlabel = "h [-]", fontfamily = "Computer Modern", legend=:bottomleft)
savefig(vmh, joinpath(HPM.BASE_DIR,"results", join(["vmh_bus_norms_","$bus_id",".pdf"])))



vmh = hcat(index, vm_ae', vm_me', vm_mm', vm_ksb', vhlim[2:end])
header = ["h" "ae" "me" "mm" "ksb" "lim"]
writedlm(joinpath(HPM.BASE_DIR,"results", join(["vmh_","$bus_id",".csv"])),  [header ; vmh], ',')

########################################
################### WITH RMS ###########
########################################
########################################
########################################

##### Processing output #######
filename = joinpath(HPM.BASE_DIR, "results", join([case, harmonic_range, "maximum efficiency.json"]))
r_me_ = Dict{String, Any}()
open(filename) do f
dicttxt = read(f, String)
global r_me_ = JSON.parse(dicttxt)
end

filename = joinpath(HPM.BASE_DIR, "results", join([case, harmonic_range, "absolute equality.json"]))
r_ae_ = Dict{String, Any}()
open(filename) do f
dicttxt = read(f, String)
global r_ae_ = JSON.parse(dicttxt)
end

filename = joinpath(HPM.BASE_DIR, "results", join([case, harmonic_range, "maximin.json"]))
r_mm_ = Dict{String, Any}()
open(filename) do f
dicttxt = read(f, String)
global r_mm_ = JSON.parse(dicttxt)
end

filename = joinpath(HPM.BASE_DIR, "results", join([case, harmonic_range, "Kalai-Smorodinsky bargaining.json"]))
r_ksb_ = Dict{String, Any}()
open(filename) do f
dicttxt = read(f, String)
global r_ksb_ = JSON.parse(dicttxt)
end



##### Harmonic current injections per node:


load_gen_buses = []

filter_gens = false

ih_me = zeros(50, length(r_me["solution"]["nw"]["2"]["load"]))
ih_ae = zeros(50, length(r_me["solution"]["nw"]["2"]["load"]))
ih_mm = zeros(50, length(r_me["solution"]["nw"]["2"]["load"]))
# ih_ksb = zeros(50, length(r_ksb["solution"]["nw"]["2"]["load"]))


for (l, load) in data["load"]
    for (g, gen) in data["gen"]
        if gen["gen_bus"] == load["load_bus"]
            push!(load_gen_buses,  gen["gen_bus"])
        end
    end
end

for (n, nw) in r_me_["solution"]["nw"]
    if n ≠ "1"
        for (l, load) in nw["load"]
            if !any(data["load"][l]["load_bus"] .== load_gen_buses) || filter_gens == false
                l_ = parse(Int, l)
                n_ = parse(Int, n)
                ih_me[n_, l_] = load["cm"]
            end
        end
    end
end

for (n, nw) in r_ae_["solution"]["nw"]
    if n ≠ "1"
        for (l, load) in nw["load"]
            if !any(data["load"][l]["load_bus"] .== load_gen_buses) || filter_gens == false
                l_ = parse(Int, l)
                n_ = parse(Int, n)
                ih_ae[n_, l_] = load["cm"]
            end
        end
    end
end

for (n, nw) in r_mm_["solution"]["nw"]
    if n ≠ "1"
        for (l, load) in nw["load"]
            if !any(data["load"][l]["load_bus"] .== load_gen_buses) || filter_gens == false
                l_ = parse(Int, l)
                n_ = parse(Int, n)
                ih_mm[n_, l_] = load["cm"]
            end
        end
    end
end

# for (n, nw) in r_ksb_["solution"]["nw"]
#     if n ≠ "1"
#         for (l, load) in nw["load"]
#             if !any(data["load"][l]["load_bus"] .== load_gen_buses) || filter_gens == false
#                 l_ = parse(Int, l)
#                 n_ = parse(Int, n)
#                 ih_ksb[n_, l_] = load["cm"]
#             end
#         end
#     end
# end

hhc_me = sum(sum(ih_me))
hhc_ae = sum(sum(ih_ae))
hhc_mm = sum(sum(ih_mm))
# hhc_ksb = sum(sum(ih_ksb))

println("HHC maximum efficiency: ", hhc_me)
println("HHC maximin: ", hhc_mm)
println("HHC absolute equality: ", hhc_ae)
# println("HHC bargaining: ", hhc_ksb)

pih_me = heatmap(log10.(ih_me)', c=cgrad(:jet, rev = true), zlabel = "\$I_{h} in pu \$", ylabel = "Load ID", xlabel = "h [-]", fontfamily = "Computer Modern")
pih_mm = heatmap(log10.(ih_mm)', c=cgrad(:jet, rev = true), zlabel = "\$I_{h} in pu \$", ylabel = "Load ID", xlabel = "h [-]", fontfamily = "Computer Modern")
pih_ae = heatmap(log10.(ih_ae)', c=cgrad(:jet, rev = true), zlabel = "\$I_{h} in pu \$", ylabel = "Load ID", xlabel = "h [-]", fontfamily = "Computer Modern")
# pih_ksb = heatmap(log10.(ih_ksb)', c=cgrad(:jet, rev = true), zlabel = "\$I_{h} in pu \$", ylabel = "Load ID", xlabel = "h [-]", fontfamily = "Computer Modern")



# zticks = ([-1 -2 -3 -4 -5 -6], ["\$10^{-1}\$","\$10^{-1}\$","\$10^{-1}\$","\$10^{-1}\$","\$10^{-1}\$","\$10^{-1}\$"])
a = rand(1000, 50)
plot(1:50,1:1000, a, linetype = :heatmap, c=cgrad(:jet, rev = true))

#outfile_ksb = joinpath(HPM.BASE_DIR,"results", join(["ih_ksb.txt"]))
outfile_ae = joinpath(HPM.BASE_DIR,"results", join(["ih_ae.txt"]))
#outfile_mm = joinpath(HPM.BASE_DIR,"results", join(["ih_mm.txt"]))
outfile_me = joinpath(HPM.BASE_DIR,"results", join(["ih_me.txt"]))

#f_ksb = open(outfile_ksb, "w")
f_ae = open(outfile_ae, "w")
f_mm = open(outfile_mm, "w")
f_me = open(outfile_me, "w")
for h in 2 : 50
    for u in 1 : 40 : 1000
 #       println(f_ksb, h, " ", u, " ", ih_ksb[h,u])
        println(f_ae, h, " ", u, " ", ih_ae[h,u])
        println(f_mm, h, " ", u, " ", ih_mm[h,u])
        println(f_me, h, " ", u, " ", ih_me[h,u])
    end
   # println(f_ksb)
    println(f_ae)
    #println(f_mm)
    println(f_me)
end
#close(f_ksb)
close(f_ae)
close(f_mm)
close(f_me)

if filter_gens == true
    savefig(pih_me, joinpath(HPM.BASE_DIR,"results", join(["ih_heat_me_filtered.pdf"])))
    savefig(pih_mm, joinpath(HPM.BASE_DIR,"results", join(["ih_heat_mm_filtered.pdf"])))
    savefig(pih_ae, joinpath(HPM.BASE_DIR,"results", join(["ih_heat_ae_filtered.pdf"])))
    savefig(pih_ksb, joinpath(HPM.BASE_DIR,"results", join(["ih_heat_ksb_filtered.pdf"])))
else
    savefig(pih_me, joinpath(HPM.BASE_DIR,"results", join(["ih_heat_me.pdf"])))
    savefig(pih_mm, joinpath(HPM.BASE_DIR,"results", join(["ih_heat_mm.pdf"])))
    savefig(pih_ae, joinpath(HPM.BASE_DIR,"results", join(["ih_heat_ae.pdf"])))
    savefig(pih_ksb, joinpath(HPM.BASE_DIR,"results", join(["ih_heat_ksb.pdf"])))
end


# zz = zeros(1, 50)
# for i in 1:50
#     zz[i] =  sqrt(Zh["$bus_id"][i]["re"]^2 + Zh["$bus_id"][i]["im"]^2)
# end
# pzh = plot(1:50, abs.(zz'), yaxis= :log10, ylim = (0.01, 100), legend = false, linewidth = 2,
# xlabel="harmonic h [-]",
# ylabel = "\$ Z_{h}\$ in pu", fontfamily = "Computer Modern", xticks = ([1,10,20,30,40,50],["1" "10" "20" "30" "40" "50"])) #, xticklabel = ["$i" for i in 1:50])
# savefig(pzh, joinpath(HPM.BASE_DIR,"results", join(["zh_bus_","$bus_id",".pdf"])))
  
    

#     ##########
# zz = zeros(length(big_buses), 50)
# for i in 1:50
#     for b in 1:length(big_buses)
#         nb = big_buses[b]
#         zz[b, i] =  sqrt(Zh["$nb"][i]["re"]^2 + Zh["$nb"][i]["im"]^2)
#     end
# end

# zh = plot(1:50, zz', yaxis= :log10, ylim = (0.01, 10), xlabel="harmonic h [-]", ylabel = "\$ Z_{h}\$ in pu", fontfamily = "Computer Modern", xticks = ([1,10,20,30,40,50],["1" "10" "20" "30" "40" "50"]))   #, xticklabel = ["$i" for i in 1:50])


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
