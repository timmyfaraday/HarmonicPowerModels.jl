# load pkgs
using Ipopt
using JuMP
using PowerModels
using HarmonicPowerModels
using SignalDecomposition
using Dierckx
using Test
using Gurobi

# pkg const
const HPM = HarmonicPowerModels
const PMs = PowerModels
const SDC = SignalDecomposition

# test functions
⪅(a,b) = (a <= b) || isapprox(a, b, atol=1e-6)

# solver
solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
solver_soc   = JuMP.optimizer_with_attributes(Gurobi.Optimizer)

# silence the warnings of PowerModels
PMs.silence()

@testset "HarmonicPowerModels.jl" begin
    
    include("hhc.jl")
    # include("hpf.jl")
    # include("hopf.jl")

end
