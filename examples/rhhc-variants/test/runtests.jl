using Test
include(joinpath(@__DIR__, "..", "compare.jl"))

@testset "Corner adversary" begin
    s = study("cases/ieee30-tom.toml")
    fixed = ParameterBox(
        s.box.names, s.box.nominal, s.box.nominal, s.box.nominal, s.box.scale)
    prepared = deepcopy(s.hd)
    nominal = solve_rhhc(prepared, s.model, fixed, s.conic;
        pce_params = deepcopy(s.pce), hpf_optimizer = s.ipopt,
        nbackground = 0, use_seeds = false)
    solution = nominal["solution"]
    result = corner_adversary(prepared, s.model, s.box, solution)
    independently_checked = [posterior_check(prepared, s.model,
                                 ParameterBox(
                                     s.box.names, point, point, point, s.box.scale), solution; n = 1)
                             for point in result.corners]
    worst = argmax([p.max_violation for p in independently_checked])
    @test result.max_violation ≈ independently_checked[worst].max_violation rtol = 1e-12
    @test result.constraint == independently_checked[worst].constraint
    @test result.count > 1
    @test any(!iszero, result.sensitivities.signs)
    stationary = corner_adversary(prepared, s.model, fixed, solution)
    direct = posterior_check(prepared, s.model, fixed, solution; n = 1)
    @test stationary.count == 1
    @test stationary.point == nominal_point(fixed)
    @test stationary.max_violation ≈ direct.max_violation rtol = 1e-12
    doubled = corner_adversary(prepared, s.model, fixed, Dict(h => 2current
    for (h, current) in solution))
    @test doubled.max_violation ≈ 4(stationary.max_violation + 1) - 1 rtol = 1e-12
end

include("stars.jl")

include("transmission.jl")

include("distribution.jl")

include("infeasible.jl")
