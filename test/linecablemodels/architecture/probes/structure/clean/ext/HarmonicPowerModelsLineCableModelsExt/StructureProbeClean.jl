module StructureProbeClean
import HarmonicPowerModels: branch_data, ParameterBox
import LineCableModels: AbstractFormulation
struct Item end
Base.length(::Item) = 1
branch_data(::AbstractFormulation, ::Val{:StructureProbeClean}) = 1
Base.getindex(::ParameterBox, ::AbstractFormulation, ::Val{:ArchitectureOwnershipProbe}) = 1
quantity(::Item) = 1
end
