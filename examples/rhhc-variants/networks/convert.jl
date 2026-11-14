using PowerImpedance, PowerModels
using LinearAlgebra: Diagonal, norm

# Run once in the source package environment. Studies read the saved JSON files.
const SOURCE = only(ARGS)
const ORIGIN = "c0547d03670b2ca4432459274c54b879d345d07b"
readchomp(`git -C $SOURCE rev-parse HEAD`) == ORIGIN || error("source commit differs")
isempty(readchomp(`git -C $SOURCE diff $ORIGIN -- examples/IEEE9_bus_impedance.jl examples/IEEE39bus.jl`)) ||
    error("source examples differ from their commit")

function source_example(file, stop)
    text = read(joinpath(SOURCE, "examples", file), String)
    prefix = first(split(text, stop))
    # Plotting is outside the network conversion.
    return replace(prefix, "using CairoMakie" => "")
end

function save_network(network, name, script, source_lines, routes, loads)
    data, bases, mapping = convert(network, ACDCPowerFlow().backend)
    normalize(data)
    data["name"] = name
    data["principle"] = "absolute equality" # `examples/shhc/data/3node_test.m:17`
    data["origin"] = Dict("commit" => ORIGIN, "script" => "examples/" * script,
        "source_lines" => source_lines,
        "conversion" => "networks/convert.jl, save_network and normalize",
        "defaults" => "PowerImpedance src/Network/powermodels_data.jl:34-53,71-97",
        "representation" => "PowerImpedance src/Network/NetworkBuilder/powerflow.jl:154-184",
        "frequency_Hz" => 50.0)
    for (key, bus) in data["bus"]
        # Convert each local phase RMS voltage base to line-to-line kV.
        # The source impedance base already includes the three-phase factor.
        # `PowerModelsACDC src/core/data.jl:15-31` defines those source bases.
        bus["base_kv"] *= sqrt(3)
        bus["standard"] = "IEC61000-3-6:2008" # `examples/shhc/data/3node_test.m:33-35`
        bus["ref_angle"] = 0.0 # same study angle convention
        bus["origin"] = "PowerImpedance examples/" * script * ":" * source_lines
    end
    for table in ("branch", "gen", "shunt"), row in values(data[table])

        row["origin"] = "PowerImpedance examples/" * script * ":" * source_lines
    end
    for gen in values(data["gen"])
        bus = data["bus"][string(gen["gen_bus"])]
        bus["bus_type"] = 4 # HPM ideal-source convention, `src/core/data.jl:59-86`
        bus["standard"] = "Clean Bus"
        # Express the source's zero scalar cost as a PowerModels polynomial.
        # `examples/shhc/data/3node_test.m:46-49` uses this zero-cost record.
        gen["model"] = 2
        gen["ncost"] = 3
        gen["cost"] = [0.0, 0.0, 0.0]
        gen["inf"] = 1
        gen["gsc"] = 0.0
        gen["bsc"] = 0.0
    end
    for (element, length_m) in routes
        id = string(mapping[element].compkey)
        data["branch"][id]["length_m"] = length_m
        data["branch"][id]["origin"] = "PowerImpedance examples/" * script * ":" *
                                       source_lines * " " *
                                       string(element)
    end
    # PowerImpedance represents grounded impedances as shunts. For HPM, retain
    # their 50 Hz power as loads so they can host the selected harmonic units.
    # Pure capacitors stay as shunts. HPM owns subsequent harmonic scaling.
    load_indices = Dict{Symbol, Int}()
    for element in loads
        id = string(mapping[element].compkey)
        shunt = pop!(data["shunt"], id)
        n = length(data["load"]) + 1
        load_indices[element] = n
        bus = shunt["shunt_bus"]
        data["load"][string(n)] = Dict{String, Any}("index" => n, "load_bus" => bus,
            "source_id" => Any["bus", bus], "status" => shunt["status"],
            "pd" => shunt["gs"], "qd" => -shunt["bs"],
            "origin" =>
                "PowerImpedance examples/" * script * ":" * source_lines * " " *
                string(element))
    end
    data["conversion_buses"] = Dict(string(row.net) => row.bus
    for row in network.topology.connections if row.bus != 0)
    # A vector records the source names without creating a PowerModels component table.
    data["conversion_elements"] = [merge(Dict("name" => string(k)),
                                       haskey(load_indices, k) ?
                                       Dict("table" => "load", "index" => load_indices[k]) :
                                       Dict("table" => v.pmtype, "index" => v.compkey))
                                   for (k, v) in pairs(mapping)]
    if name == "ieee39"
        # MATPOWER 8.1, `data/case39.m:136`: Pmax in MW, Qmax and Qmin in MVAr.
        # https://github.com/MATPOWER/matpower/blob/8.1/data/case39.m#L136
        limits = (Pmax = 1100.0, Qmax = 300.0, Qmin = -100.0)
        rating = hypot(limits.Pmax, max(abs(limits.Qmax), abs(limits.Qmin))) /
                 data["baseMVA"]
        row = data["branch"][string(mapping[:Zg39].compkey)]
        # At nominal voltage, the per-unit current rating equals the power rating.
        for field in ("rate_a", "rate_b", "rate_c", "c_rating_a")
            row[field] = rating
        end
        row["rating_origin"] = "MATPOWER tag 8.1, data/case39.m:136, " *
                               "https://github.com/MATPOWER/matpower/blob/8.1/data/case39.m#L136; " *
                               "Pmax=1100 MW, Qmax=300 MVAr, Qmin=-100 MVAr; " *
                               "S=hypot(Pmax,max(abs(Qmax),abs(Qmin))) MVA; " *
                               "rate_a=rate_b=rate_c=c_rating_a=S/baseMVA at nominal voltage"
    end
    destination = joinpath(@__DIR__, name * ".json")
    PowerModels.export_file(destination, data)
    restored = PowerModels.parse_file(destination; validate = false)
    for table in ("bus", "branch", "gen", "shunt", "load")
        isequal(data[table], restored[table]) || error("saved $name $table differs")
    end
    println(name, ": ", length(data["bus"]), " buses, ", length(data["branch"]),
        " branches, ", length(data["load"]), " loads")
end

# PowerImpedance sets nominal voltage levels through symmetric 0.9/1.1 limits:
# `src/Network/Components/transformer/transformer.jl:365-369` and
# `src/Network/Components/source/source.jl:184-186` in the source checkout.
# Let k be the local nominal voltage divided by the source global voltage base.
# V_old = diag(k) V_new, so Y_new = diag(k) Y_old diag(k).
function normalize(data)
    factors = Dict(id => (bus["vmin"] + bus["vmax"]) / 2
    for (id, bus) in data["bus"])
    all(==(1.0), values(factors)) && return data
    before = PowerModels.calc_admittance_matrix(data)
    for (id, bus) in data["bus"]
        k = factors[id]
        bus["base_kv"] *= k
        for field in ("vm", "vmin", "vmax")
            bus[field] /= k
        end
    end
    for gen in values(data["gen"])
        gen["vg"] /= factors[string(gen["gen_bus"])]
    end
    for row in values(data["shunt"])
        k = factors[string(row["shunt_bus"])]
        row["gs"] *= k^2
        row["bs"] *= k^2
    end
    for row in values(data["branch"])
        kf, kt = factors[string(row["f_bus"])], factors[string(row["t_bus"])]
        row["br_r"] /= kt^2
        row["br_x"] /= kt^2
        for field in ("g_fr", "b_fr", "g_to", "b_to")
            row[field] *= kt^2
        end
        row["tap"] *= kt / kf
    end
    after = PowerModels.calc_admittance_matrix(data)
    diagonal = Diagonal([factors[string(id)] for id in before.idx_to_bus])
    expected = diagonal * before.matrix * diagonal
    isapprox(after.matrix, expected; rtol = 1e-12) ||
        error("local voltage bases changed the physical admittance")
    println("Local voltage-base admittance identity: relative error ",
        norm(after.matrix - expected) / norm(expected))
    return data
end

# IEEE9_bus_impedance.jl:18-97 supplies the passive network and its connections.
include_string(Main, source_example("IEEE9_bus_impedance.jl", "network = NetworkBuilder.define"))
network9 = NetworkBuilder.define(elements, connections;
    options = (; voltageBase = 230 / sqrt(3))) # source IEEE9 line 175 voltage base
save_network(network9, "ieee9", "IEEE9_bus_impedance.jl", "18-97",
    # Original routes at lines 31-48, in meters.
    (:tl78 => 76e3, :tl89 => 106e3, :tl75 => 170e3,
        :tl96 => 179e3, :tl54 => 89e3, :tl64 => 97e3), (:l5, :l6, :l8))

# IEEE39bus.jl:18-187 supplies passive elements before the STATCOM study.
include_string(Main, source_example("IEEE39bus.jl", "# ## The STATCOM"))
# Reuse the source's connection rows, omitting the separate STATCOM ports.
source39 = read(joinpath(SOURCE, "examples", "IEEE39bus.jl"), String)
rows39 = "generator_rows =" * split(split(source39, "generator_rows ="; limit = 2)[2],
    "connections = (statcom_connections"; limit = 2)[1]
include_string(Main, rows39)
network39 = NetworkBuilder.define(
    elements, (generator_rows..., branch_rows..., load_rows...);
    options = (; voltageBase = 345 / sqrt(3))) # source IEEE39bus.jl:335-336
save_network(network39, "ieee39", "IEEE39bus.jl", "18-187,297-324",
    (row.name => row.length for row in line_data),
    (row.name for row in load_data if row.P > 0))
