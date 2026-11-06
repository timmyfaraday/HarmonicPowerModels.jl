include(joinpath(@__DIR__, "builders.jl"))

# Olex, printed pp. 12 and 13. The table is in `data/register/README.md`, catalogue values used as checks.
# Olex resistance uses a conductor temperature of 90 °C.
# The reactance and capacitance refer to 50 Hz.
const OLEX = Dict("U33-400" => (0.101, 0.116, 0.267),
    "U33-500" => (0.0797, 0.111, 0.297), "U33-630" => (0.0629, 0.108, 0.329),
    "U33-800" => (0.0507, 0.104, 0.366))

function checkmodel(s, specname, formname, temperature)
    spec = specification(specname)
    family_name = spec.family_name
    family = get(s.families, family_name, (;))
    bindings = (;
        (k => Symbol(v) for (k, v) in pairs(get(s.config.bindings, family_name, (;))))...)
    fn = problem(
        spec, family, formation(formname), s.mat, s.config.earth, temperature, bindings)
    model = ExactPiModel(Dict(1 => (problem = fn, len = 1.0, zbase = 1.0)),
        CoaxialFormulation(transposition = true), s.box; f1 = s.config.f1)
    return model.branches[1].zy
end

function checks(io, s)
    p = nominal_point(s.box)
    println(io, "Data checks. DC uses the positive-sequence real impedance at 0.01 Hz and 20 °C.")
    println(io, "All values below are per km. R and X are in Ω/km. C is in µF/km.")
    println(io,
        "spec | DC LCM | catalogue R20 | R50 at 20 °C | R50 at 90 °C | Olex R90 | X50 | Olex X | C50 | Olex C | diameter mm | catalogue mm")
    for name in ("TOM-630", "U33-400", "U33-500", "U33-630", "U33-800")
        spec = specification(name)
        form = name == "TOM-630" ? "trefoil-tom" : "trefoil"
        c = checkmodel(s, name, form, 20.0)
        low = c(p, 0.01).positive
        v = c(p, 50.0).positive
        family = name == "TOM-630" ? s.families.TOM : s.families.U33
        bindings = name == "TOM-630" ? get(s.config.bindings, :TOM, (;)) :
                   get(s.config.bindings, :U33, (;))
        q = coordinates(family, (; (k => Symbol(v) for (k, v) in pairs(bindings))...), p)
        outer = 2r_ex(design(spec.family, spec, q, s.mat).outer) * 1e3
        if haskey(OLEX, name)
            hot = checkmodel(s, name, form, 90.0)(p, 50.0).positive
            ref = OLEX[name]
            @printf(io,
                "%s | %.6f | %.6f | %.6f | %.6f | %.6f | %.6f | %.6f | %.6f | %.6f | %.3f | %.3f\n",
                name, 1e3real(low.z), 1e3spec.R20, 1e3real(v.z), 1e3real(hot.z), ref[1],
                1e3imag(v.z), ref[2], 1e9imag(v.y) / (2π * 50), ref[3], outer,
                1e3spec.diameter)
            spacing = outer / 1e3
            screen_diameter = spec.d_c + 2(spec.t_sc1 + spec.t_ins + spec.t_sc2) +
                              spec.d_wire
            omega = 2π * 50
            screen_resistance = spec.R20_screen * (1 + s.mat.copper.alpha * (90 - 20))
            # IEC 60287-1-1:2006, Section 2.3.1 gives the mutual reactance and circulating-current loss factor.
            mutual = omega * 2e-7 * log(2spacing / screen_diameter)
            loss = screen_resistance * mutual^2 / (screen_resistance^2 + mutual^2)
            corrected_r = ref[1] + 1e3loss
            classical_x = omega * (0.05 + 0.2log(2spacing / spec.d_c)) * 1e-3
            reduction = 1e3mutual^3 / (screen_resistance^2 + mutual^2)
            corrected_x = classical_x - reduction
            @printf(io,
                "U33 bonding check %s: R90 Olex+screen %.6f, LCM %.6f, difference %+.3f %%; X50 classical %.6f, minus screen %.6f, LCM %.6f, difference %+.3f %% [Ω/km].\n",
                name, corrected_r, 1e3real(hot.z), 100(1e3real(hot.z) / corrected_r - 1),
                classical_x, corrected_x, 1e3imag(v.z), 100(1e3imag(v.z) / corrected_x - 1))
        else
            @printf(io,
                "%s | %.6f | unavailable | %.6f | unavailable | unavailable | %.6f | unavailable | %.6f | unavailable | %.3f | unavailable\n",
                name, 1e3real(low.z), 1e3real(v.z), 1e3imag(v.z),
                1e9imag(v.y) / (2π * 50), outer)
        end
    end
    println(io,
        "ACSR: spec | DC LCM | catalogue R20 | R50 at 20 °C | X50 LCM | source X50 | C50 LCM | scenario C50")
    for name in ("Panther", "Lion", "Bear", "Goat")
        spec = specification(name)
        c = checkmodel(s, name, "O132-D", 20.0)
        v = c(p, 50.0).positive
        @printf(io, "%s | %.6f | %.6f | %.6f | %.6f | %.6f | %.6f | %.6f\n",
            name, 1e3real(c(p, 0.01).positive.z), 1e3specification(name).R20,
            1e3real(v.z), 1e3imag(v.z), OHL_X[132.0], 1e9imag(v.y) / (2π * 50), 0.0)
        # BS EN 50182:2001, Section 5.10 excludes galvanized steel from nominal DC resistance.
        steel = spec.n_steel * π * (spec.d_steel / 2)^2 / s.mat.s1a.rho
        aluminum = spec.n_aluminium * π * (spec.d_aluminium / 2)^2 / s.mat.al1.rho
        @printf(io,
            "ACSR steel conductance share %s: %.6f %% (wire areas and resistivities).\n",
            name, 100steel / (steel + aluminum))
    end
    gw = checkmodel(s, "GW-7x3.15", "O132-D", 20.0)
    @printf(io,
        "GW-7x3.15 DC LCM %.6f Ω/km. Register O132.gw.R20: 3.55 Ω/km, interval 3.30–4.01 Ω/km.\n",
        1e3real(gw(p, 0.01).positive.z))
    println(io, "The earth-wire check uses three copies as phases at the same mean-height formation.")
    println(io, "At 0.01 Hz their real positive-sequence impedance approaches the wire's DC resistance.")
end

function scan(path)
    s = study(path)
    name = splitext(basename(s.file))[1]
    mkpath(joinpath(STUDIES, "output"))
    dest = joinpath(STUDIES, "output", name * "-sequences.csv")
    p = nominal_point(s.box)
    # HPM modal_scan defaults: hmin 1, step 0.02, hmax 1.25 times the highest harmonic.
    hs = sort(unique([collect(1.0:0.02:(1.25maximum(s.config.harmonics)));
                      collect(s.config.harmonics)]))
    ids = circuits(s)
    open(dest, "w") do io
        println(io,
            "circuit,f_Hz,h,z0_real_ohm_per_m,z0_imag_ohm_per_m,z1_real_ohm_per_m,z1_imag_ohm_per_m,y0_real_S_per_m,y0_imag_S_per_m,y1_real_S_per_m,y1_imag_S_per_m")
        for b in ids, h in hs

            c = s.model.branches[b].zy
            v = c(p, h * s.config.f1)
            println(io,
                join(
                    (s.details[b].circuit, h * s.config.f1, h,
                        real(v.zero.z), imag(v.zero.z), real(v.positive.z), imag(v.positive.z),
                        real(v.zero.y), imag(v.zero.y), real(v.positive.y), imag(v.positive.y)),
                    ','))
        end
    end
    report = joinpath(STUDIES, "output", name * "-scan.txt")
    open(report, "w") do io
        header(io, s)
        println(io, "Sequence file: ", basename(dest), "; ",
            length(hs), " frequencies per circuit.")
        println(io, "branch,buses,kV,spec,formation,length_m,circuit,source_X_ohm_per_km,scenario_C_uF_per_km")
        for b in sort(collect(keys(s.details)))
            d = s.details[b]
            println(io,
                join(
                    (b, "$(d.f_bus)–$(d.t_bus)", d.kv, d.spec, d.formation,
                        d.len, d.circuit, d.source_x, 1e6d.source_c),
                    ','))
        end
        println(io, "50 Hz sequences per meter: circuit,z0,z1,y0,y1")
        for b in ids
            v = s.model.branches[b].zy(p, s.config.f1)
            println(io, join(
                (
                    s.details[b].circuit, v.zero.z, v.positive.z, v.zero.y, v.positive.y), ','))
        end
        resonances = modal_scan(deepcopy(s.hd), s.model, p)
        println(io, "resonances at the nominal point: ",
            join([@sprintf("h = %.4f", r.h) for r in resonances], ", "))
        flush(io)
        checks(io, s)
    end
    print(read(report, String))
end

length(ARGS) == 1 || error("usage: scan.jl cases/<case>.toml")
scan(only(ARGS))
