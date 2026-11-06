"""
    Circuit{B, F, P, V}

Serve a three-phase LCM circuit per sequence. Calls at the box's point type and a
`Float64` frequency share stored results. Other point or frequency types compute directly.
"""
struct Circuit{B <: Function, F <: AbstractFormulation, P <: NamedTuple, V <: NamedTuple}
    "User function from a point and frequency \\[Hz\\] to an LCM problem."
    problem::B
    "LCM formulation built once and passed unchanged to each computation."
    formulation::F
    "Three-phase sequence transformation."
    fortescue::Fortescue
    "Sequences in \\[Ω/m\\] and \\[S/m\\], keyed by point and frequency \\[Hz\\]."
    sequences::Dict{Tuple{P, Float64}, V}
    "Lock that guards access to the sequences."
    lock::ReentrantLock

    @doc """
        Circuit(problem::Function, formulation::AbstractFormulation, box::ParameterBox)

    Construct a circuit whose `problem(p, f)` returns an LCM `LineParametersProblem`
    at frequency `f` \\[Hz\\]. The stored point type comes from `box`.
    """
    function Circuit(problem::Function, formulation::AbstractFormulation, box::ParameterBox)
        P = typeof(nominal_point(box))
        V = NamedTuple{(:zero, :positive, :negative),
            Tuple{
                NamedTuple{(:z, :y), Tuple{ComplexF64, ComplexF64}},
                NamedTuple{(:z, :y), Tuple{ComplexF64, ComplexF64}},
                NamedTuple{(:z, :y), Tuple{ComplexF64, ComplexF64}}}}
        return new{typeof(problem), typeof(formulation), P, V}(
            problem, formulation, Fortescue(3), Dict{Tuple{P, Float64}, V}(), ReentrantLock())
    end
end

"""
    compute(c::Circuit, p::NamedTuple, f::Real)

Compute LCM's matrices at point `p` and frequency `f` \\[Hz\\], transform and
validate them, and return their diagonal as `zero`, `positive` and `negative`.
Each member contains series impedance `z` \\[Ω/m\\] and shunt admittance `y` \\[S/m\\].
"""
function compute(c::Circuit, p::NamedTuple, f::Real)
    result = compute(c.problem(p, f), c.formulation)
    Z = c.fortescue(result.Z[:, :, 1])
    Y = c.fortescue(result.Y[:, :, 1])
    validate(Z, c)
    validate(Y, c)
    return (zero = (z = Z[1, 1], y = Y[1, 1]),
        positive = (z = Z[2, 2], y = Y[2, 2]),
        negative = (z = Z[3, 3], y = Y[3, 3]))
end

"""
    validate(M::AbstractMatrix, c::Circuit)

Return `M` unchanged. Warn when its largest off-diagonal magnitude exceeds
`c.fortescue.tol` times its largest diagonal magnitude, with a floor of machine
precision. This is LCM's `v0.0.0:src/utils/Utils.jl` `offdiag_ratio` check.
"""
function validate(M::AbstractMatrix, c::Circuit)
    d = zero(real(float(eltype(M))))
    o = zero(d)
    for j in axes(M, 2), i in axes(M, 1)

        if i == j
            d = max(d, abs(M[i, j]))
        else
            o = max(o, abs(M[i, j]))
        end
    end
    if o > c.fortescue.tol * max(d, eps(real(float(eltype(M)))))
        @warn "Sequence matrix is not diagonal within tolerance" ratio=o/max(d, eps(real(float(eltype(M))))) tol=c.fortescue.tol
    end
    return M
end

"""
    (c::Circuit{B, F, P})(p::P, f::Float64) where {B<:Function, F<:AbstractFormulation, P<:NamedTuple}

Return all sequences at point `p` and frequency `f` \\[Hz\\]. Store the first
written result under a lock. Compute outside the lock when the key is absent.
The members `z` and `y` have units \\[Ω/m\\] and \\[S/m\\].
"""
function (c::Circuit{B, F, P})(p::P, f::Float64) where {B <: Function, F <: AbstractFormulation, P <: NamedTuple}
    @lock c.lock begin
        haskey(c.sequences, (p, f)) && return c.sequences[(p, f)]
    end
    result = compute(c, p, f)
    @lock c.lock begin
        return get!(c.sequences, (p, f), result)
    end
end

"""
    (c::Circuit)(p::NamedTuple, f::Real)

Compute all sequences directly. Preserve the types of `p` and
frequency `f` \\[Hz\\]. Return `z` \\[Ω/m\\] and `y` \\[S/m\\] in each sequence.
"""
(c::Circuit)(p::NamedTuple, f::Real) = compute(c, p, f)

"""
    (c::Circuit)(p, f, ::ZeroSequence)

Return zero-sequence `(z, y)` in \\[Ω/m\\] and \\[S/m\\] at `p` and `f` \\[Hz\\].
"""
(c::Circuit)(p, f, ::ZeroSequence) = c(p, f).zero

"""
    (c::Circuit)(p, f, ::PositiveSequence)

Return positive-sequence `(z, y)` in \\[Ω/m\\] and \\[S/m\\] at `p` and `f` \\[Hz\\].
"""
(c::Circuit)(p, f, ::PositiveSequence) = c(p, f).positive

"""
    (c::Circuit)(p, f, ::NegativeSequence)

Return negative-sequence `(z, y)` in \\[Ω/m\\] and \\[S/m\\] at `p` and `f` \\[Hz\\].
"""
(c::Circuit)(p, f, ::NegativeSequence) = c(p, f).negative
