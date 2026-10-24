################################################################################
# HarmonicPowerModels.jl — capacitive scenarios for the IEEE30 recombination    #
# study (feature/shhc-recombination, Phase 1)                                    #
#                                                                               #
# Each scenario modifies the parsed fundamental data of                          #
# data/30bus_network_shhc_pv.m in place, before `replicate`:                     #
#   base     — the Phase 0 baseline, no capacitor banks                          #
#   banks    — the original IEEE30 bus shunts restored: 19 MVAr at bus 10 and    #
#              4.3 MVAr at bus 24                                                #
#   cable33  — all 33 kV lines converted to underground cable                    #
#   cable132 — all 132 kV lines converted to underground cable                   #
#   cableall — all lines converted to underground cable                          #
#   cable132c, cableallc — as above, with the cable charging compensated by     #
#              shunt reactors at both cable ends (K_COMP of the susceptance),    #
#              as cable systems are operated in practice; uncompensated, the    #
#              132 kV cables overload branch 6-8 at the fundamental              #
#                                                                               #
# Cable conversion. The route length of each line is recovered from its         #
# reactance assuming an overhead line, and the line is replaced by an XLPE      #
# cable of the same length. Typical per-km values (assumptions, not data):      #
#   132 kV  overhead x' = 0.40 Ω/km;  cable r' = 0.03 Ω/km, x' = 0.13 Ω/km,     #
#           c' = 0.20 µF/km                                                     #
#    33 kV  overhead x' = 0.35 Ω/km;  cable r' = 0.10 Ω/km, x' = 0.11 Ω/km,     #
#           c' = 0.30 µF/km                                                     #
# The transformer branches are left unchanged. The total cable susceptance      #
# ωc'L is split equally over both ends of the π-model.                          #
#                                                                               #
# Operating point. Every scenario is re-dispatched with the AC-OPF recipe of    #
# the baseline — bus 1 at 1.0 p.u., slack reactive limits lifted, branch        #
# ratings enforced, lower voltage limit relaxed to 0.90 p.u. — with the upper   #
# voltage limit tightened to 1.05 p.u. (OPF only) so the operating point keeps  #
# RMS voltage headroom for the harmonics; the PV setpoints Pg and Vg are taken  #
# from it. For "base" this reproduces the setpoints in the data file (it never  #
# reaches 1.05 p.u.). Without the re-dispatch the cable scenarios overload      #
# branch 1-2 or 6-8 at the fundamental; without the 1.05 p.u. cap they sit at   #
# 1.06 p.u. and have no RMS headroom.                                           #
# Author: Tom Van Acker                                                          #
################################################################################

const F_NOM = 50.0

const OHL_X = Dict(132.0 => 0.40, 33.0 => 0.35)                      # Ω/km
const K_COMP = 0.8                                                   # reactor compensation
const CABLE = Dict(132.0 => (r = 0.03, x = 0.13, c = 0.20e-6),       # Ω/km, Ω/km, F/km
                    33.0 => (r = 0.10, x = 0.11, c = 0.30e-6))

"Restore the original IEEE30 bus shunts (capacitor banks) at buses 10 and 24."
function add_original_banks!(data)
    data["shunt"] = Dict{String,Any}(
        "1" => Dict{String,Any}("index" => 1, "shunt_bus" => 10, "gs" => 0.0,
                                "bs" => 0.19,  "status" => 1, "source_id" => ["bus", 10]),
        "2" => Dict{String,Any}("index" => 2, "shunt_bus" => 24, "gs" => 0.0,
                                "bs" => 0.043, "status" => 1, "source_id" => ["bus", 24]))
    return data
end

"""
Convert every non-transformer line at the voltage levels `kvs` to underground
cable. With `comp > 0`, a shunt reactor of `comp` times the cable's half
susceptance is placed at each cable end (inductive bus shunt, scaled by 1/h at
the harmonics).
"""
function to_cable!(data, kvs; comp = 0.0)
    sbase = data["baseMVA"]
    ω     = 2π * F_NOM
    for br in values(data["branch"])
        br["transformer"] && continue
        kv = data["bus"]["$(br["f_bus"])"]["base_kv"]
        kv in kvs || continue
        zb  = kv^2 / sbase                         # Ω
        len = br["br_x"] * zb / OHL_X[kv]          # km, from the overhead reactance
        cab = CABLE[kv]
        br["br_r"] = cab.r * len / zb
        br["br_x"] = cab.x * len / zb
        b          = ω * cab.c * len * zb          # p.u., total line charging
        br["b_fr"] = b / 2
        br["b_to"] = b / 2
        br["cable_km"] = len
        if comp > 0
            shunts = get!(data, "shunt", Dict{String,Any}())
            for bus in (br["f_bus"], br["t_bus"])
                k = length(shunts) + 1
                shunts["$k"] = Dict{String,Any}("index" => k, "shunt_bus" => bus, "gs" => 0.0,
                                                "bs" => -comp * b / 2, "status" => 1,
                                                "source_id" => ["reactor", k])
            end
        end
    end
    return data
end

"""
Remove the charging of every line (non-transformer branch): a purely inductive
network, as the inductive case of Van Acker and Ergun (PESGM 2025). Capacitive
elements are then only those a scenario adds explicitly.
"""
function inductive!(data)
    for br in values(data["branch"])
        br["transformer"] && continue
        br["b_fr"] = 0.0; br["b_to"] = 0.0
    end
    return data
end

"Add one capacitor bank of `bs` p.u. at `bus`."
function add_bank!(data, bus::Int, bs::Real)
    shunts = get!(data, "shunt", Dict{String,Any}())
    k = length(shunts) + 1
    shunts["$k"] = Dict{String,Any}("index" => k, "shunt_bus" => bus, "gs" => 0.0, "bs" => bs,
                                    "status" => 1, "source_id" => ["bus", bus])
    return data
end

# paper scenarios: purely inductive base network, then capacitive elements added
const PAPER_BANK_BUS  = 29          # bus of the capacitor bank (positive impact, see sens_ieee30.jl)
const PAPER_BANK_SIZE = 0.19        # p.u. (19 Mvar)

const SCENARIOS = Dict{String,Function}(
    "base"     => d -> d,
    "banks"    => add_original_banks!,
    "cable33"  => d -> to_cable!(d, (33.0,)),
    "cable132" => d -> to_cable!(d, (132.0,)),
    "cableall" => d -> to_cable!(d, (33.0, 132.0)),
    "cable132c" => d -> to_cable!(d, (132.0,); comp = K_COMP),
    "cableallc" => d -> to_cable!(d, (33.0, 132.0); comp = K_COMP),
    "inductive"           => d -> inductive!(d),
    "inductive_bank"      => d -> add_bank!(inductive!(d), PAPER_BANK_BUS, PAPER_BANK_SIZE),
    "inductive_cable33"   => d -> to_cable!(inductive!(d), (33.0,)),
    "inductive_cable132c" => d -> to_cable!(inductive!(d), (132.0,); comp = K_COMP),
)

"""
Re-dispatch `data` with the baseline AC-OPF recipe and write the resulting
active power and voltage setpoints into the generators (`pg`, `vg`).
"""
function redispatch!(data, nlp_optimizer; vmin = 0.90, lift_q = :slack)
    d = deepcopy(data)
    for (i, bus) in d["bus"]
        bus["bus_type"] == 3 ? (bus["vmin"] = bus["vmax"] = 1.0) :
                               (bus["vmin"] = vmin; bus["vmax"] = min(bus["vmax"], 1.05))
    end
    for gen in values(d["gen"])
        (lift_q == :all || d["bus"]["$(gen["gen_bus"])"]["bus_type"] == 3) &&
            (gen["qmin"] = -10.0; gen["qmax"] = 10.0)
    end
    res = PowerModels.solve_ac_opf(d, nlp_optimizer)
    st  = res["termination_status"]
    st in (PowerModels.LOCALLY_SOLVED, PowerModels.OPTIMAL) ||
        error("redispatch!: AC-OPF did not solve ($st)")
    for (g, gen) in data["gen"]
        gen["pg"] = res["solution"]["gen"][g]["pg"]
        gen["vg"] = res["solution"]["bus"]["$(gen["gen_bus"])"]["vm"]
    end
    return data
end

"Parse the single-slack IEEE30 case, apply scenario `name` and re-dispatch it."
function ieee30_scenario(name::String, nlp_optimizer)
    haskey(SCENARIOS, name) || error("unknown scenario \"$name\"; use one of $(sort(collect(keys(SCENARIOS))))")
    data = HarmonicPowerModels.parse_file(joinpath(@__DIR__, "data", "30bus_network_shhc_pv.m"))
    SCENARIOS[name](data)
    # paper scenarios (purely inductive base): original voltage limits 0.94-1.05 p.u.,
    # with the reactive limits of all generators and synchronous condensers lifted at
    # the fundamental, since the line charging no longer supports the voltage
    startswith(name, "inductive") ? redispatch!(data, nlp_optimizer; vmin = 0.94, lift_q = :all) :
                                    redispatch!(data, nlp_optimizer)
    return data
end

"Parse the single-slack IEEE30 case with one capacitor bank `bs` [p.u.] at `bus`, re-dispatched (legacy, with line charging)."
function ieee30_bank(bus::Int, bs::Real, nlp_optimizer)
    data = HarmonicPowerModels.parse_file(joinpath(@__DIR__, "data", "30bus_network_shhc_pv.m"))
    data["shunt"] = Dict{String,Any}("1" => Dict{String,Any}("index" => 1, "shunt_bus" => bus, "gs" => 0.0,
                                     "bs" => bs, "status" => 1, "source_id" => ["bus", bus]))
    redispatch!(data, nlp_optimizer)
    return data
end

"Purely inductive IEEE30 case with one capacitor bank `bs` [p.u.] at `bus`, re-dispatched as the paper scenarios."
function ieee30_inductive_bank(bus::Int, bs::Real, nlp_optimizer)
    data = HarmonicPowerModels.parse_file(joinpath(@__DIR__, "data", "30bus_network_shhc_pv.m"))
    add_bank!(inductive!(data), bus, bs)
    redispatch!(data, nlp_optimizer; vmin = 0.94, lift_q = :all)
    return data
end
