using HarmonicPowerModels, LineCableModels, PowerModels, JuMP, Ipopt, Clarabel
using Printf, LinearAlgebra, TOML

const HPM = HarmonicPowerModels
const STUDIES = @__DIR__
include(joinpath(STUDIES, "..", "shhc", "pce_params.jl"))
include(joinpath(STUDIES, "..", "shhc", "ieee30_scenarios.jl"))
PowerModels.silence()
BLAS.set_num_threads(1)

# Decode each TOML record once. The files keep the quantities and their citations.
record(x) = x
record(x::AbstractVector) = Tuple(record(v) for v in x)
function record(x::AbstractDict)
    haskey(x, "value") && return record(x["value"])
    return (; (Symbol(k) => record(v) for (k, v) in x)...)
end
readrecord(path) = record(TOML.parsefile(path))

function materials(path)
    rows = readrecord(path)
    return (; (k => merge(v, (kind = Val(Symbol(v.kind)),)) for (k, v) in pairs(rows))...)
end

function metal(row; rho = row.rho)
    return Material(:conductor, rho, get(row, :eps_r, 1.0), get(row, :mu_r, 1.0),
        get(row, :T0, 20.0), row.alpha)
end
function dielectric(row; eps_r = row.eps_r)
    Material(:insulator, row.rho, eps_r;
        tan_delta = get(row, :tan_delta, 0.0))
end
semiconductor(row) = Material(:semicon, row.rho, row.eps_r)

function specification(name)
    row = readrecord(joinpath(STUDIES, "data", "specs", name * ".toml"))
    return merge(row, (family = Val(Symbol(row.family)), family_name = Symbol(row.family)))
end

function formation(name)
    row = readrecord(joinpath(STUDIES, "data", "formations", name * ".toml"))
    return merge(row,
        (name = name, kind = Val(Symbol(row.kind)),
            reference = Val(Symbol(get(row, :reference, "support")))))
end

function coordinates(family, bindings, p)
    return (;
        (k => get(p, get(bindings, k, k), row.nominal)
    for (k, row) in pairs(family))...)
end

function design(::Val{:TOM}, spec, q, mat)
    aluminum = metal(mat.tom_aluminium)
    copper = metal(mat.tom_copper)
    semicon = semiconductor(mat.tom_semicon)
    xlpe = dielectric(mat.tom_xlpe; eps_r = q.eps_xlpe)
    pe = dielectric(mat.tom_pe)
    return @design "TOM-630" begin
        @terminal :core begin
            core(aluminum; r = q.r_core)
        end
        screen(semicon; t = spec.t_sc1, tag = :semicon_inner)
        insulation(xlpe; t = q.t_xlpe)
        screen(semicon; t = spec.t_sc2, tag = :semicon_outer)
        @terminal :sheath begin
            sheath(copper; t = spec.t_tube)
        end
        jacket(pe; t = spec.t_jkt)
    end
end

function design(::Val{:U33}, spec, q, mat)
    radius = (spec.d_c + q.dd) / 2
    # GUSTAVSEN2005 III-A: R20 and diameter jointly define the equivalent solid core.
    aluminum = metal(mat.aluminium; rho = q.kR * spec.R20 * π * radius^2)
    copper = metal(mat.copper)
    inner = semiconductor(mat.semicon_inner)
    outer = semiconductor(mat.semicon_outer)
    xlpe = dielectric(mat.xlpe; eps_r = q.eps)
    pvc = dielectric(mat.pvc)
    screen_radius = radius + spec.t_sc1 + q.kt * spec.t_ins + spec.t_sc2
    # TB531 4.2.2.5 and GUSTAVSEN2016 conclusion 4: tube area matches screen R20.
    tube = sqrt(screen_radius^2 + mat.copper.rho / (π * spec.R20_screen)) - screen_radius
    thickness = q.kj * spec.t_jkt + spec.d_wire - tube
    return @design "U33" begin
        @terminal :core begin
            core(aluminum; r = radius)
        end
        screen(inner; t = spec.t_sc1, tag = :semicon_inner)
        insulation(xlpe; t = q.kt * spec.t_ins)
        screen(outer; t = spec.t_sc2, tag = :semicon_outer)
        @terminal :sheath begin
            sheath(copper; t = tube)
        end
        jacket(pvc; t = thickness)
    end
end

function design(::Val{:ACSR}, spec, q, mat)
    # The family factor acts on the terminal's resistance through both metals' resistivity.
    steel = metal(mat.s1a; rho = q.kR * mat.s1a.rho)
    aluminum = metal(mat.al1; rho = q.kR * mat.al1.rho)
    return @design "ACSR" begin
        @terminal :core begin
            stranded(steel; center = Circle(spec.d_steel / 2),
                shape = Circle(spec.d_steel / 2), n = (spec.n_steel - 1,),
                lay = LayRatio(spec.lay_steel), boundary = Circle(3spec.d_steel / 2),
                tag = :steel)
            wires(aluminum; shape = Circle(spec.d_aluminium / 2), n = spec.courses[1],
                lay = LayRatio(spec.lay_inner), tag = :aluminium_inner)
            wires(aluminum; shape = Circle(spec.d_aluminium / 2), n = spec.courses[2],
                lay = LayRatio(spec.lay_outer), tag = :aluminium_outer)
        end
    end
end

function design(::Val{:steel}, spec, q, mat)
    steel = metal(mat.s1a)
    return @design "GW-7x3.15" begin
        @terminal :core begin
            stranded(steel; center = Circle(spec.diameter / 2),
                shape = Circle(spec.diameter / 2), n = (spec.count - 1,),
                lay = LayRatio(spec.lay), boundary = Circle(3spec.diameter / 2))
        end
    end
end

centre(::Val{:centre}, row, spacing) = Placement(0.0, -row.depth)
centre(::Val{:lower}, row, spacing) = Placement(0.0, -row.depth + spacing / (2sqrt(3)))

function system(::Val{:trefoil}, conductor, row, mat)
    spacing = 2r_ex(conductor.outer)
    positions = trefoil(conductor; spacing, at = centre(row.reference, row, spacing))
    return build(CableSystem, fill(conductor, 3), positions;
        connections = [(core = k, sheath = 0) for k in 1:3], id = "trefoil", line_length = 1.0)
end

function system(::Val{:overhead}, conductor, row, mat)
    # The parabolic span's mean height is support height minus two thirds of the sag.
    positions = [Placement(x, h - (2 / 3) * row.phase_sag)
                 for (x, h) in zip(row.phase_x, row.phase_height)]
    if !haskey(row, :earth_wire)
        return build(CableSystem, fill(conductor, 3), positions;
            connections = [(core = k,) for k in 1:3], id = row.name, line_length = 1.0)
    end
    push!(positions, Placement(row.shield_x, row.shield_height - (2 / 3) * row.shield_sag))
    return build(CableSystem, [conductor, conductor, conductor, row.shield], positions;
        connections = [(core = 1,), (core = 2,), (core = 3,), (core = 0,)],
        id = row.name, line_length = 1.0)
end

arrangement(::Val{:trefoil}, row, mat) = row
function arrangement(::Val{:overhead}, row, mat)
    haskey(row, :earth_wire) || return row
    shield_spec = specification(row.earth_wire)
    return merge(row, (shield = design(shield_spec.family, shield_spec, (;), mat),))
end

function problem(spec, family, row, mat, earth, temperature, bindings)
    layout = arrangement(row.kind, row, mat)
    soil_binding = Symbol(earth.binding)
    # Each circuit gives its branches the same problem function.
    return (p, f) -> begin
        q = coordinates(family, bindings, p)
        conductor = design(spec.family, spec, q, mat)
        circuit = system(layout.kind, conductor, layout, mat)
        soil = homogeneous(rho = get(p, soil_binding, earth.rho), eps_r = earth.eps_r)
        LineParametersProblem(circuit, soil; temperature, frequencies = [f])
    end
end

function parameterbox(rows, families, names)
    entries = map(names) do text
        name = Symbol(text)
        row = getproperty(rows, name)
        values = haskey(row, :family) ?
                 getproperty(getproperty(families, Symbol(row.family)),
            Symbol(row.coordinate)) : row
        return name =>
            (values.nominal, values.lower, values.upper, Symbol(get(row, :scale, "linear")))
    end
    return ParameterBox(entries...)
end

function isxfmr(b)
    get(b, "transformer", false) || !(get(b, "tap", 1.0) in (0.0, 1.0)) ||
        get(b, "shift", 0.0) != 0.0
end

function transformers(rows, data)
    for (name, row) in pairs(rows)
        id = string(name)
        haskey(data["branch"], id) || error("declared transformer branch $id is absent")
        br = data["branch"][id]
        isxfmr(br) || error("declared branch $id is not a transformer")
        (br["f_bus"], br["t_bus"]) == (row.f_bus, row.t_bus) ||
            error("transformer $id does not join declared buses $(row.f_bus)–$(row.t_bus)")
    end
    return rows
end

# MATPOWER route comments remain the source of lengths, in meters after decoding.
function network(config, ipopt)
    source = config.network
    source.kind == "scenario" && return ieee30_scenario(source.scenario, ipopt)
    file = normpath(joinpath(STUDIES, source.file))
    data = HPM.parse_file(file)
    if source.kind == "matpower"
        text = read(file, String)
        table = match(r"mpc\.branch\s*=\s*\[(.*?)\];"s, text)
        for (id, line) in enumerate(filter(
            l -> !isempty(strip(l)) && !startswith(strip(l), "%"),
            split(table.captures[1], '\n')))
            route = match(r"(?:[lL]\s*=\s*|[,;]\s*)([0-9.]+)\s*km", line)
            if route !== nothing
                data["branch"][string(id)]["length_m"] = 1e3parse(Float64, route.captures[1])
            end
        end
    end
    return data
end

function route_length(br, kv, zbase, config)
    haskey(br, "length_m") && return br["length_m"]
    haskey(br, "cable_km") && return 1e3br["cable_km"]
    config.network.kind == "scenario" && return 1e3br["br_x"] * zbase / OHL_X[kv]
    error("LCM branch $(br["index"]) has no route length in its source data")
end

function harmonic_units(count::Integer, data)
    sort(parse.(Int, collect(keys(data["load"]))))[1:count]
end
harmonic_units(ids::Tuple, data) = collect(ids)

function study(path)
    file = isfile(path) ? abspath(path) : joinpath(STUDIES, path)
    config = readrecord(file)
    mat = materials(joinpath(STUDIES, "data", "materials.toml"))
    families = (;
        (Symbol(name) => readrecord(joinpath(STUDIES, "data", "families", name * ".toml"))
    for name in ("TOM", "U33", "ACSR"))...)
    box = parameterbox(config.box, families, config.box_names)
    ipopt = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
    conic = optimizer_with_attributes(Clarabel.Optimizer, "verbose" => false)
    data = network(config, ipopt)
    haskey(config, :principle) && (data["principle"] = config.principle)
    transformers(get(config, :transformers, (;)), data)
    definitions = Dict{Tuple{String, String}, Function}()
    details = Dict{Int, NamedTuple}()
    branches = Dict{Int, NamedTuple}()
    for (name, row) in pairs(get(config, :branches, (;)))
        b = parse(Int, string(name))
        br = data["branch"][string(b)]
        kv = data["bus"][string(br["f_bus"])]["base_kv"]
        zbase = kv^2 / data["baseMVA"]
        len = route_length(br, kv, zbase, config)
        key = (row.spec, row.formation)
        if !haskey(definitions, key)
            spec = specification(row.spec)
            family_name = spec.family_name
            definitions[key] = problem(spec, getproperty(families, family_name),
                formation(row.formation), mat, config.earth, config.temperature,
                (;
                    (k => Symbol(v)
                for (k, v) in pairs(get(config.bindings, family_name, (;))))...))
        end
        branches[b] = (problem = definitions[key], len = len, zbase = zbase)
        details[b] = (branch = b, f_bus = br["f_bus"], t_bus = br["t_bus"], kv = kv,
            spec = row.spec, formation = row.formation, len = len,
            circuit = row.spec * "/" * row.formation,
            source_x = br["br_x"] * zbase / (len / 1e3),
            source_c = (br["b_fr"] + br["b_to"]) / (2π * config.f1 * zbase * (len / 1e3)))
    end
    if config.model == "exact"
        formulation = CoaxialFormulation(transposition = true)
        model = ExactPiModel(branches, formulation, box; f1 = config.f1)
        p = nominal_point(box)
        for b in keys(branches), (k, v) in pairs(branch_data(model, b, p, 1.0, :positive))

            data["branch"][string(b)][string(k)] = v
        end
    end
    units = harmonic_units(config.units, data)
    all(id -> haskey(data["load"], string(id)), units) ||
        error("harmonic unit is absent from the source load table")
    hd = HPM.replicate(deepcopy(data); H = [1; collect(config.harmonics)])
    if config.model == "scaled"
        classes = Dict(b => Symbol(name) for (name, ids) in pairs(config.classes)
        for b in ids)
        model = ScaledBranchModel(hd, classes)
    end
    hd["angle_reference"] = "local"
    for h in config.harmonics, k in collect(keys(hd["nw"][string(h)]["load"]))

        parse(Int, k) in units || delete!(hd["nw"][string(h)]["load"], k)
    end
    pce = Dict{String, Any}(deepcopy(PCE_PARAMS)...)
    pce["harmonics"] = Dict{Int, Any}(h => (PCE_PARAMS["harmonics"][h].v_min,
                                          PCE_PARAMS["harmonics"][h].v_max)
    for h in config.harmonics)
    pce["harmonic_units"] = units
    pce["bound"] = config.bound
    pce["epsilon"] = config.epsilon
    return (; file, config, mat, families, box, model, details, hd, pce, ipopt, conic)
end

function circuits(s)
    ids = sort(collect(keys(s.details)))
    return unique(b -> s.details[b].circuit, ids)
end

function header(io, s)
    println(io, "git commit: ", readchomp(`git -C $STUDIES rev-parse HEAD`))
    println(io, "case: ", relpath(s.file, dirname(dirname(STUDIES))))
    println(io, "threads: ", Threads.nthreads(), "; solver: Clarabel")
end

function declarations(io, s)
    isempty(get(s.config, :transformers, (;))) && return
    println(io, "Transformer declarations. Not applied to the network: HPM reads vector groups")
    println(io, "from xfmr, which has no tap. IEEE30 transformers are tapped branches.")
    for (id, row) in sort(collect(pairs(s.config.transformers)); by = x -> parse(Int, string(first(x))))
        br = s.hd["nw"]["1"]["branch"][string(id)]
        kv1 = s.hd["nw"]["1"]["bus"][string(br["f_bus"])]["base_kv"]
        kv2 = s.hd["nw"]["1"]["bus"][string(br["t_bus"])]["base_kv"]
        zb1, zb2 = kv1^2 / s.hd["nw"]["1"]["baseMVA"], kv2^2 / s.hd["nw"]["1"]["baseMVA"]
        println(io, "branch ", id, " (", row.f_bus, "–", row.t_bus, "): vg=", row.vg,
            " gnd1=", row.gnd1, " gnd2=", row.gnd2,
            " re1=", row.re1_ohm / zb1, " xe1=", row.xe1_ohm / zb1,
            " re2=", row.re2_ohm / zb2, " xe2=", row.xe2_ohm / zb2, " [pu]")
        println(io, row.note)
        haskey(row, :resistance_note) && println(io, row.resistance_note)
    end
end
