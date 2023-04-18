# load pkgs
using Ipopt
using PowerModels
using HarmonicPowerModels
using SignalDecomposition
using Dierckx
using Test

# pkg const
const HPM = HarmonicPowerModels
const PMs = PowerModels
const SDC = SignalDecomposition

# solver
solver = Ipopt.Optimizer

@testset "HarmonicPowerModels.jl" begin
    
    include("hpf.jl")
    include("hopf.jl")
    include("xfmr.jl")
end
