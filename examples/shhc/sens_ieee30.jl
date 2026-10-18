################################################################################
# HarmonicPowerModels.jl — IEEE30 sensitivities for the SHHC paper              #
#                                                                               #
#   1. bank location: one 19 Mvar bank (0.19 p.u.) at each bus in turn,         #
#      re-dispatched; populations A and AB, N = 21                               #
#   2. bank size at bus 10: 0 … 0.40 p.u.; populations A and AB, N = 21          #
#   3. number of converter types T = 1 … 5, mean firing delays 0, 15, 30, 45,   #
#      60° assigned cyclically, one shared germ per controlled type (30° wide);  #
#      base and banks, N = 21                                                    #
# Reported: exact 95% hosting (binding IHD), its equivalent summation exponent  #
# per harmonic, and the robust bound.                                            #
#                                                                               #
# Output: output/ieee30_sens_location.csv, _size.csv, _types.csv                #
# Author: Tom Van Acker                                                          #
################################################################################

include(joinpath(@__DIR__, "recomb_ieee30.jl"))
include(joinpath(@__DIR__, "converter_types.jl"))

const NS = 21
out(name) = open(joinpath(@__DIR__, "output", "ieee30_sens_$(name).csv"), "w")
row(io, pre, r) = (println(io, pre, ",", join([@sprintf("%.6e", sum(r[h][k] for h in harm)) for k in (:q95, :rob, :det, :arith)], ","),
                           ",", join([string(r[h][:q95_alpha]) for h in harm], ",")); flush(io))
header = "q95,rob,det,arith,alpha5,alpha7,alpha11,alpha13"

io = out("location"); println(io, "bus,population,", header)
for bus in 2:30, pop in ("A", "AB")
    try
        data  = ieee30_bank(bus, 0.19, ipopt)
        units = sort(parse.(Int, collect(keys(data["load"]))))[1:NS]
        r = recombination("bank$bus", NS, "local"; data = data, types = pop == "AB" ? typeB_germ(units) : nothing)
        row(io, "$bus,$pop", r)
    catch err
        println(io, "$bus,$pop,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN"); flush(io)
    end
end
close(io)

io = out("size"); println(io, "bs,population,", header)
for bs in (0.0, 0.05, 0.10, 0.15, 0.19, 0.25, 0.30, 0.40), pop in ("A", "AB")
    try
        data  = ieee30_bank(10, bs, ipopt)
        units = sort(parse.(Int, collect(keys(data["load"]))))[1:NS]
        r = recombination("size$bs", NS, "local"; data = data, types = pop == "AB" ? typeB_germ(units) : nothing)
        row(io, "$bs,$pop", r)
    catch err
        println(io, "$bs,$pop,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN"); flush(io)
    end
end
close(io)

io = out("types"); println(io, "scenario,T,", header)
for scen in ("base", "banks"), T in 1:5
    data  = ieee30_scenario(scen, ipopt)
    units = sort(parse.(Int, collect(keys(data["load"]))))[1:NS]
    delays = [0.0, 15.0, 30.0, 45.0, 60.0][1:T]
    r = recombination(scen, NS, "local"; data = data, types = T == 1 ? nothing : types_cyclic(units, delays))
    row(io, "$scen,$T", r)
end
close(io)
println("✓ ieee30_sens_*.csv")
