@testset "Fortescue" begin
    F = Fortescue(3)
    @test F.matrix * F.matrix' ≈ I
    a = cis(2π / 3)
    A = [1 1 1; 1 a^2 a; 1 a a^2]
    M = rand(ComplexF64, 3, 3)
    M = (M + transpose(M)) / 2
    @test F(M) ≈ inv(A) * M * A

    ms, mm = 2.0 + 3.0im, 0.3 + 0.7im
    M = fill(mm, 3, 3) + (ms - mm) * I
    result = F(M)
    @test diag(result) ≈ [ms + 2mm, ms - mm, ms - mm]
    @test maximum(abs, result - Diagonal(diag(result))) < 1e-14

    F4 = Fortescue(4)
    @test F4.matrix * F4.matrix' ≈ I
    terms = [2.0 + 3.0im, 0.2 + 0.4im, 0.1 + 0.3im, 0.2 + 0.4im]
    M4 = [terms[mod(j - i, 4) + 1] for i in 1:4, j in 1:4]
    result4 = F4(M4)
    @test diag(result4) ≈ [sum(terms), terms[1] - terms[3],
        terms[1] - 2terms[2] + terms[3], terms[1] - terms[3]]
    @test maximum(abs, result4 - Diagonal(diag(result4))) < 1e-14
    @test_throws DimensionMismatch("matrix size (6, 6) differs from Fortescue size (3, 3)") F(zeros(6, 6))
end
