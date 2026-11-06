using HarmonicPowerModels: HarmonicPowerModels, ExactPiModel, ParameterBox,
                           nominal_point, branch_data, latent
using LineCableModels: LineCableModels, Material, @design, @terminal, core,
                       screen, insulation, sheath, jacket, r_ex, trefoil, flat, Placement,
                       build,
                       CableSystem, homogeneous, LineParametersProblem, CoaxialFormulation,
                       compute
using LinearAlgebra: I, diag, Diagonal, BLAS
using Test: @testset, @test, @test_throws, @test_logs, @test_nowarn, @inferred

const Ext = Base.get_extension(HarmonicPowerModels, :HarmonicPowerModelsLineCableModelsExt)
using .Ext: Circuit, Fortescue

BLAS.set_num_threads(1)
include("fixtures.jl")

@testset "LCM extension" begin
    include("sequence.jl")
    include("circuit.jl")
    include("protocol.jl")
end
