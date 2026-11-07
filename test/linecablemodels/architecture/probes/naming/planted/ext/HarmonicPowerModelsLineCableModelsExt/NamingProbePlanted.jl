module NamingProbePlanted
struct Circuit
    matrix::Matrix{Float64}
    Circuit(matrix::Matrix{Float64}) = new(matrix)
end
unknown_name(x) = x
line_matrix(x) = x + 1
public line_matrix
end
