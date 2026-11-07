module NamingProbeClean
struct Circuit
    matrix::Matrix{Float64}
    Circuit(matrix::Matrix{Float64}) = new(matrix)
end
matrix(c::Circuit) = c.matrix
end
