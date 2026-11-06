include(joinpath(@__DIR__, "builders.jl"))

function sweep(path)
    s = study(path; all_coordinates = true)
    p = nominal_point(s.box)
    name = splitext(basename(s.file))[1]
    dest = joinpath(STUDIES, "output", name * "-sweep.txt")
    mkpath(dirname(dest))
    admitted = String[]
    open(dest, "w") do io
        header(io, s)
        println(io, "Relative changes compare each endpoint with the nominal value.")
        println(io,
            "coordinate | maximum relative change | circuit | sequence | quantity | Hz | end | decision")
        for (index, key) in enumerate(s.box.names)
            row = getproperty(s.config.box, key)
            best = (change = -1.0, circuit = "", sequence = :zero,
                quantity = :z, frequency = 50.0, endpoint = "lower")
            for endpoint in ("lower", "upper")
                value = getproperty(s.box, Symbol(endpoint))[index]
                point = merge(p, NamedTuple{(key,)}((value,)))
                for b in circuits(s)
                    family = string(specification(s.details[b].spec).family_name)
                    reads = row.family == "EARTH" || row.family == family ||
                            (haskey(formation(s.details[b].formation), :structure) &&
                             row.family ==
                             string(formation(s.details[b].formation).structure))
                    reads || continue
                    circuit = s.model.branches[b].zy
                    for frequency in (50.0, 650.0)
                        nominal = circuit(p, frequency)
                        result = circuit(point, frequency)
                        for sequence in (:zero, :positive, :negative), quantity in (:z, :y)

                            reference = getproperty(getproperty(nominal, sequence), quantity)
                            changed = getproperty(getproperty(result, sequence), quantity)
                            relative = abs(changed - reference) / abs(reference)
                            if relative > best.change
                                best = (change = relative, circuit = s.details[b].circuit,
                                    sequence, quantity, frequency, endpoint)
                            end
                        end
                    end
                end
            end
            decision = best.change >= 1e-3 ? "yes" : "no"
            decision == "yes" && push!(admitted, string(key))
            @printf(io, "%s | %.12g | %s | %s | %s | %.0f | %s | %s\n",
                key, best.change, best.circuit, best.sequence, best.quantity,
                best.frequency, best.endpoint, decision)
            flush(io)
        end
        println(io, "box_names = [", join(repr.(admitted), ", "), "]")
        spec = specification("Bear")
        q = coordinates(s.families.ACSR, (;), (;))
        corner = merge(q,
            (d_aluminium = s.families.ACSR.d_aluminium.upper,
                d_steel = s.families.ACSR.d_steel.upper))
        realized = realize(spec, corner, q)
        conductor = design(spec.family, spec, corner, s.mat; realized)
        diameter = 2r_ex(conductor.outer)
        @assert spec.diameter.lower <= diameter <= spec.diameter.upper
        @assert spec.mass.lower <= realized.mass <= spec.mass.upper
        @assert realize(spec, q, q).factor == 1
        println(io, "Bear projection check: upper aluminum and steel wire diameters, nominal lay ratios.")
        @printf(io,
            "t = %.12g; diameter = %.12g m (%+.9f percent); mass = %.12g kg/m (%+.9f percent).\n",
            realized.factor, diameter, 100(diameter / spec.diameter.nominal - 1),
            realized.mass, 100(realized.mass / spec.mass.nominal - 1))
        @printf(io, "Effective wires: aluminum %.12g m, steel %.12g m.\n", realized.da,
            realized.ds)
        println(io, "EN 50182: diameter within +/-1 percent and mass within +/-2 percent; passed.")
        println(io, "Inner lay upper ", s.families.ACSR.lay_inner.upper,
            " equals steel lay lower ", s.families.ACSR.lay_steel.lower,
            "; that order cannot bind except at equality.")
        for code in ("Bear", "Goat", "Lion", "Panther"),
            a in (s.families.ACSR.d_aluminium.lower, s.families.ACSR.d_aluminium.upper),
            d in (s.families.ACSR.d_steel.lower, s.families.ACSR.d_steel.upper),
            inner in (s.families.ACSR.lay_inner.lower, s.families.ACSR.lay_inner.upper),
            outer in (s.families.ACSR.lay_outer.lower, s.families.ACSR.lay_outer.upper),
            steel in (s.families.ACSR.lay_steel.lower, s.families.ACSR.lay_steel.upper)
            product = specification(code)
            point = merge(q,
                (d_aluminium = a, d_steel = d,
                    lay_inner = inner, lay_outer = outer, lay_steel = steel))
            result = realize(product, point, q)
            @assert product.diameter.lower <= result.diameter <= product.diameter.upper
            @assert product.mass.lower <= result.mass <= product.mass.upper
            @assert result.outer <= result.inner <= steel
        end
        println(io,
            "All 128 wire-diameter and lay endpoint combinations satisfy both aggregate limits and lay order.")
    end
    print(read(dest, String))
end

length(ARGS) == 1 || error("usage: sweep.jl cases/<case>.toml")
sweep(only(ARGS))
