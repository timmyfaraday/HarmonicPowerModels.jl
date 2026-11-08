@testset "Antwerp physical cables" begin
    s = study("cases/antwerp-physical.toml")
    source = study("cases/antwerp.toml")
    p = nominal_point(s.box)
    @test length(s.details) == 7
    @test length(s.lcm) == 1
    # `sweep.jl cases/antwerp-physical.toml`, November 8, 2026, HPM 6059a62, LCM e180b435e.
    @test length(s.box) == 12
    @test s.hd["nw"]["1"]["xfmr"] == source.hd["nw"]["1"]["xfmr"]
    @test count(t -> t["vg"] == "Yy0", values(s.hd["nw"]["1"]["xfmr"])) == 4
    @test count(t -> t["vg"] == "Yd11", values(s.hd["nw"]["1"]["xfmr"])) == 7
    for b in (1, 2, 10, 11, 12, 13, 14)
        @test s.hd["nw"]["1"]["branch"][string(b)] ==
              source.hd["nw"]["1"]["branch"][string(b)]
    end
    circuit = s.model.branches[3].zy
    spec = specification("U155-Cu1600")
    q = coordinates(s.families.U155, (;), (;))
    cable = design(spec.family, spec, q, s.mat)
    # Candidate study section 3.5: U155.diameter and U155.screen_count.
    @test 2r_ex(cable.outer) ≈ 0.109 rtol = 1e-12
    @test spec.screen_count == 60
    cold = merge(p, (theta_U155 = 20.0,))
    @test real(circuit(cold, 0.01).positive.z) ≈ q.R20 rtol = 1e-4
    # medium_network_shhc.m:101-107, individual source routes and ratings.
    for (b, len) in zip(3:9, (10600.0, 10600.0, 1917.0, 1810.0, 1850.0, 2039.0, 2039.0))
        @test s.details[b].len == len
        @test s.model.branches[b].zy === circuit
        @test s.hd["nw"]["1"]["branch"][string(b)]["rate_a"] == 3.0
        for f in (50.0, 250.0, 350.0, 550.0, 650.0)
            v = circuit(p, f)
            @test v.negative.z ≈ v.positive.z rtol = 1e-12
            @test v.negative.y ≈ v.positive.y rtol = 1e-12
            @test imag(v.positive.y) > 0
            pi = branch_data(s.model, b, p, f / 50, :positive)
            gamma = sqrt(v.positive.z * v.positive.y)
            characteristic = sqrt(v.positive.z / v.positive.y)
            expected = characteristic * sinh(gamma * len) / s.model.branches[b].zbase
            @test complex(pi.br_r, pi.br_x) ≈ expected rtol = 1e-12
            @test pi.b_fr ≈
                  imag(tanh(gamma * len / 2) / characteristic) *
                  s.model.branches[b].zbase rtol = 1e-12
            @test pi.b_fr == pi.b_to
        end
    end
    # Reference run: `builders.jl`, `study("cases/antwerp-physical.toml")`, November 8, 2026.
    # HPM parent 6059a62, LCM e180b435e. Circuit at 50 Hz and the two temperature endpoints.
    # Total AC resistance includes circulating screen losses and falls slightly at these endpoints.
    @test real(circuit(cold, 50.0).positive.z) ≈ 3.0014191694832257e-5 rtol = 1e-12
    @test real(circuit(p, 50.0).positive.z) ≈ 2.994829746284598e-5 rtol = 1e-12
    @test real(circuit(cold, 50.0).positive.z) > real(circuit(p, 50.0).positive.z)
    fixed = ParameterBox(
        s.box.names, s.box.nominal, s.box.nominal, s.box.nominal, s.box.scale)
    prepared = deepcopy(s.hd)
    nominal = solve_rhhc(prepared, s.model, fixed, s.conic;
        pce_params = deepcopy(s.pce), hpf_optimizer = s.ipopt, nbackground = 0, use_seeds = false)
    result = solve_rhhc(prepared, s.model, s.box, s.conic;
        pce_params = deepcopy(s.pce), hpf_optimizer = s.ipopt, posterior = s.config.posterior)
    @test nominal["signed_off"]
    @test posterior_check(prepared, s.model, fixed, nominal["solution"]; n = 1).max_violation <=
          0
    @test result["signed_off"]
    @test result["posterior"].violated_fraction == 0
    @test result["posterior"].max_violation <= 0
    @test all(br -> max(br["cm_fr"], br["cm_to"]) <= br["c_rating"],
        values(prepared["nw"]["1"]["branch"]))
end

@testset "IEEE30 cable sequence references" begin
    # `scan.jl cases/ieee30-full.toml`, November 8, 2026, HPM 6059a62, LCM e180b435e.
    s = study("cases/ieee30-full.toml")
    p = nominal_point(s.box)
    checked = Set{String}()
    for line in eachline(joinpath(STUDIES, "output", "ieee30-full-sequences.csv"))
        fields = split(line, ',')
        startswith(fields[1], "U33-") || continue
        f = parse(Float64, fields[2])
        f in (50.0, 250.0, 350.0, 550.0, 650.0) || continue
        b = first(b for (b, row) in s.details if row.circuit == fields[1])
        actual = s.model.branches[b].zy(p, f)
        reference = parse.(Float64, fields[4:end])
        @test isequal(actual.zero.z, complex(reference[1], reference[2]))
        @test isequal(actual.positive.z, complex(reference[3], reference[4]))
        @test isequal(actual.zero.y, complex(reference[5], reference[6]))
        @test isequal(actual.positive.y, complex(reference[7], reference[8]))
        push!(checked, fields[1])
    end
    @test length(checked) == 4
end
