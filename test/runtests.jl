# load pkgs
using Ipopt
using PowerModels
using HarmonicPowerModels
using SignalDecomposition
using Test

# pkg const
const HPM = HarmonicPowerModels
const PMs = PowerModels
const SDC = SignalDecomposition

# solver
solver = Ipopt.Optimizer

@testset "HarmonicPowerModels.jl" begin
    
    include("hpf.jl")
    include("xfmr.jl")

end
