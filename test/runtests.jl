# load pkgs
using Ipopt
using PowerModels
using HarmonicPowerModels
using Test

# pkg const
const _PMs = PowerModels
const _HPM = HarmonicPowerModels

# solver
solver = Ipopt.Optimizer

@testset "HarmonicPowerModels.jl" begin
    
    include("xfmr.jl")

end
