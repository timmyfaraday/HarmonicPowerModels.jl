################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Authors: Tom Van Acker, Frederik Geth, Hakan Ergun                           #
################################################################################
# Changelog:                                                                   #
################################################################################

# load pkgs
using Test

using JuMP
using PowerModels
using HarmonicPowerModels

using Clarabel
using Ipopt

# pkg const
const PMs = PowerModels
const HPM = HarmonicPowerModels

# test functions
≈(a,b) = isapprox(a, b, atol=1e-6)
⪅(a,b) = (a <= b) || isapprox(a, b, atol=1e-6)

# solvers
solver_nlp = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
solver_soc = JuMP.optimizer_with_attributes(Clarabel.Optimizer, "verbose" => 0)

# silence the warnings of PowerModels
PMs.silence()

@testset "HarmonicPowerModels.jl" begin
    
    # models
    include("hhc.jl")
    include("hpf.jl")
    include("hopf.jl")

end
