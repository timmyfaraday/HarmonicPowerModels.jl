using Test
using JuMP
using Clarabel
using HarmonicPowerModels

@testset "test_interval" begin
    model = Model(Clarabel.Optimizer)
    set_silent(model)
    @variable(model, -i <= x[i=1:2] <= i)
    @objective(model, Max, HarmonicPowerModels.min_box(
        model,
        [MOI.Interval(1.0, 1.0), MOI.Interval(2.0, 3.0)],
        x,
        MOI.Interval(-1.0, 0.0),
    ))
    optimize!(model)
    @test value.(x) ≈ [1, 2] rtol=1e-6
    @test objective_value(model) ≈ 4 rtol=1e-6
end
