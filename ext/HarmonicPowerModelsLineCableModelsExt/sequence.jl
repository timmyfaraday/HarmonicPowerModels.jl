"""
    AbstractSequence

Represent a sequence of HPM's protocol, as defined by `harmonic_sequence` in
HPM's `src/core/data.jl`.
"""
abstract type AbstractSequence end

"""
    ZeroSequence

Represent the zero sequence of HPM's `harmonic_sequence` (`src/core/data.jl`),
for harmonic orders divisible by three.
"""
struct ZeroSequence <: AbstractSequence end

"""
    PositiveSequence

Represent the positive sequence of HPM's `harmonic_sequence` (`src/core/data.jl`),
for harmonic orders congruent to one modulo three.
"""
struct PositiveSequence <: AbstractSequence end

"""
    NegativeSequence

Represent the negative sequence of HPM's `harmonic_sequence` (`src/core/data.jl`),
for harmonic orders congruent to two modulo three.
"""
struct NegativeSequence <: AbstractSequence end

"""
    Fortescue

Represent the generalized Fortescue transformation, ported from LCM
`v0.0.0:src/engine/transforms/fortescue.jl`. For `n` phases, its unitary matrix is

```math
F_{km} = \\frac{\\exp(\\mathrm{j}2\\pi km/n)}{\\sqrt{n}},
\\qquad k,m=0,\\ldots,n-1, \\qquad F^{-1}=F^{\\mathsf{H}}.
```

For `n = 3`, ``F M F^{\\mathsf{H}} = A^{-1} M A``, with the classical matrix
```math
A = \\begin{bmatrix} 1 & 1 & 1 \\\\ 1 & a^2 & a \\\\ 1 & a & a^2 \\end{bmatrix},
```

with
``a=\\exp(\\mathrm{j}2\\pi/3)``. The diagonal order is zero, positive, negative.
For a circulant impedance with self term `zs` and mutual term `zm`, both in
\\[Ω/m\\], ``z_0=zs+2zm`` and ``z_1=z_2=zs-zm``. The same law holds for
shunt admittance in \\[S/m\\].
"""
struct Fortescue
    "Unitary matrix of the N-point DFT."
    matrix::Matrix{ComplexF64}
    "Relative tolerance of the circuit's diagonal check."
    tol::Float64

    @doc """
        Fortescue(n::Integer; tol::Real = 1e-4)

    Construct the transformation for `n` phases with relative tolerance `tol`.
    Throw `ArgumentError` when `n < 1`.
    """
    function Fortescue(n::Integer; tol::Real = 1e-4)
        n >= 1 || throw(ArgumentError("n = $n must be at least 1"))
        return new([cis(2π * k * m / n) / sqrt(n) for k in 0:(n - 1), m in 0:(n - 1)], Float64(tol))
    end
end

"""
    (F::Fortescue)(M::AbstractMatrix)

Return ``F M_s F^{\\mathsf{H}}``, where ``M_s=(M+M^{\\mathsf{T}})/2`` enforces
reciprocity. The result retains the units of `M`, either \\[Ω/m\\] or \\[S/m\\].
Refuse a matrix whose size differs from `F.matrix`.
"""
function (F::Fortescue)(M::AbstractMatrix)
    validate(M, F)
    Ms = (M + transpose(M)) / 2
    return F.matrix * Ms * F.matrix'
end

"""
    validate(M::AbstractMatrix, F::Fortescue)

Return `M` unchanged, or throw `DimensionMismatch` with both matrix sizes.
"""
function validate(M::AbstractMatrix, F::Fortescue)
    size(M) == size(F.matrix) || throw(DimensionMismatch(
        "matrix size $(size(M)) differs from Fortescue size $(size(F.matrix))"))
    return M
end
