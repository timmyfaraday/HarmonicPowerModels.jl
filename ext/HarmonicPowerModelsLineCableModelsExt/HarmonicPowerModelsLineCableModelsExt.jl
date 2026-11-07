"""
    HarmonicPowerModelsLineCableModelsExt

Serve impedance per unit length from LCM \\[Ω/m\\] and shunt admittance \\[S/m\\] to HPM's
sequence protocol through the generalized Fortescue transformation.
"""
module HarmonicPowerModelsLineCableModelsExt

import HarmonicPowerModels: ExactPiModel, ParameterBox, nominal_point
import LineCableModels: compute, validate, select, AbstractFormulation

include("sequence.jl")
include("circuit.jl")
include("entry.jl")

public Circuit, Fortescue, AbstractSequence, ZeroSequence, PositiveSequence,
       NegativeSequence

end
