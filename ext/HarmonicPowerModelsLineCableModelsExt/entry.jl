"""
    select(::Type{AbstractSequence}, sequence::Symbol)
    select(::Type{AbstractSequence}, ::Val{S}) where {S}

Lift HPM's sequence to its type. Refuse a sequence other than `:zero`, `:positive`
or `:negative` with `ArgumentError` with the sequence.
"""
select(::Type{AbstractSequence}, sequence::Symbol) = select(AbstractSequence, Val(sequence))
select(::Type{AbstractSequence}, ::Val{:zero}) = ZeroSequence()
select(::Type{AbstractSequence}, ::Val{:positive}) = PositiveSequence()
select(::Type{AbstractSequence}, ::Val{:negative}) = NegativeSequence()
function select(::Type{AbstractSequence}, ::Val{S}) where {S}
    throw(ArgumentError("sequence :$S is unavailable; Circuit serves :zero, :positive and :negative"))
end

"""
    (c::Circuit)(p, f, sequence::Symbol)

Serve HPM's sequence at `p` and frequency `f` \\[Hz\\]. Return series
impedance `z` \\[Ω/m\\] and shunt admittance `y` \\[S/m\\].
"""
(c::Circuit)(p, f, sequence::Symbol) = c(p, f, select(AbstractSequence, sequence))

"""
    ExactPiModel(branches, formulation, box; f1 = 50.0)

Construct HPM's exact π-model through the extension's entry point. Each branch
record has `problem(p, f)`, which returns an LCM problem at frequency `f` \\[Hz\\],
`len` \\[m\\] as a number or function of `p`, and `zbase` \\[Ω\\]. Branches whose
problem functions are identical (`===`) share one `Circuit` and its sequence table.

LCM receives the formulation unchanged. The fundamental
frequency `f1` defaults to 50 \\[Hz\\]. The circuit supplies impedance per unit length
\\[Ω/m\\] and shunt admittance \\[S/m\\]. HPM forms the branch data in per unit.
"""
function ExactPiModel(branches::AbstractDict{Int, <:NamedTuple},
        formulation::AbstractFormulation, box::ParameterBox; f1::Real = 50.0)
    circuits = IdDict{Function, Circuit}()
    for br in values(branches)
        if !haskey(circuits, br.problem)
            circuits[br.problem] = Circuit(br.problem, formulation, box)
        end
    end
    return ExactPiModel(
        Dict(b => (zy = circuits[br.problem], len = br.len, zbase = br.zbase)
        for (b, br) in branches);
        f1 = f1)
end
