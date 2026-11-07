module StructureProbePlanted
import HarmonicPowerModels: branch_data, HarmonicPowerModels
import LineCableModels: AbstractFormulation, compute
struct Item end
Base.length(::Item) = 1
branch_data(::AbstractFormulation, ::Val{:StructureProbePlanted}) = 1
branch_data(::Val{:ArchitectureOwnershipProbe}) = 1
branch_data(::HarmonicPowerModels.ParameterBox, ::Val{:ArchitectureOwnershipProbe}) = 1
compute(::AbstractFormulation, ::Val{:ArchitectureOwnershipProbe}) = 1
Base.length(::AbstractFormulation, ::Val{:ArchitectureOwnershipProbe}) = 1
struct ParameterBox end
struct CoaxialFormulation end
function length end
include("../../test/linecablemodels/misplaced.jl")
end
