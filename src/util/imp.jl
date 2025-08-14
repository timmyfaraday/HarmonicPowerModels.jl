################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Authors: Tom Van Acker                                                       #
################################################################################
# Changelog:                                                                   #
# v0.3.0 - init                                                                #
################################################################################

# util #########################################################################
""
function add!(i::Vector{Int}, j::Vector{Int}, v::Vector{ComplexF64}, I::Int, J::Int, V::Complex)
    push!(i, I)
    push!(j, J)
    push!(v, V)
end
""
function find_sparse_matrix_idx(A::_SPA.SparseMatrixCSC, i::Int, j::Int)
    nzr     = _SPA.nzrange(A,j)
    rows    = view(_SPA.rowvals(A),nzr)
    return nzr[searchsortedfirst(rows,i)]
end

# edge #########################################################################
## branch ######################################################################
""
admittance_branch_series_ondiag(r::Float64, x::Float64, h::Float64) = 
    1 / (sqrt(h) * r + im * h * x)
""
admittance_branch_series_offdiag(r::Float64, x::Float64, h::Float64) =
   -1 / (sqrt(h) * r + im * h * x)
""
admittance_branch_shunt_ondiag(b::Float64, g::Float64, h::Float64) = 
   g / sqrt(h) + im * h * b

""
init_branch_data(branch::Dict{String,Any}) = 
   (Nb          = length(branch),
    idx_s_fr    = zeros(Int, length(branch)),
    idx_s_to    = zeros(Int, length(branch)),
    idx_sh_fr   = zeros(Int, length(branch)),
    idx_sh_to   = zeros(Int, length(branch)),
    fr          = Int[br["f_bus"] for br in values(branch)],
    to          = Int[br["t_bus"] for br in values(branch)],
    r           = Float64[br["br_r"] for br in values(branch)],
    x           = Float64[br["br_x"] for br in values(branch)],
    b_fr        = Float64[br["b_fr"] for br in values(branch)],
    g_fr        = Float64[br["g_fr"] for br in values(branch)],
    b_to        = Float64[br["b_to"] for br in values(branch)],
    g_to        = Float64[br["g_to"] for br in values(branch)])
""
function fill_branch_idx!(branch::NamedTuple, y::_SPA.SparseMatrixCSC)
    for nb in 1:branch[:Nb]
        branch[:idx_s_fr][nb]   = find_sparse_matrix_idx(y, branch[:fr][nb], branch[:to][nb])
        branch[:idx_s_to][nb]   = find_sparse_matrix_idx(y, branch[:to][nb], branch[:fr][nb])
        branch[:idx_sh_fr][nb]  = find_sparse_matrix_idx(y, branch[:fr][nb], branch[:fr][nb])
        branch[:idx_sh_to][nb]  = find_sparse_matrix_idx(y, branch[:to][nb], branch[:to][nb])
end end

## xfmr ########################################################################
""
admittance_xfmr_series_ondiag(r::Float64, x::Float64, h::Float64) = 
    1 / (sqrt(h) * r + im * h * x)
""
admittance_xfmr_series_offdiag(r::Float64, x::Float64, h::Float64) =
   -1 / (sqrt(h) * r + im * h * x)
""
admittance_xfmr_shunt_ondiag(b::Float64, g::Float64, h::Float64) = 
   g / sqrt(h) + im * b / h
""
init_xfmr_data(xfmr::Dict{String,Any}) = 
   (Nx          = length(xfmr),
    idx_s_fr    = zeros(Int, length(xfmr)),
    idx_s_to    = zeros(Int, length(xfmr)), 
    idx_sh_to   = zeros(Int, length(xfmr)),
    fr          = Int[xf["f_bus"] for xf in values(xfmr)], 
    to          = Int[xf["t_bus"] for xf in values(xfmr)],
    r           = Float64[xf["r1"] + xf["r2"] for xf in values(xfmr)],
    x           = Float64[xf["xsc"] for xf in values(xfmr)],
    b           = zeros(Float64, length(xfmr)),
    g           = Float64[xf["gsh"] for xf in values(xfmr)])
""
function fill_xfmr_idx!(xfmr::NamedTuple, y::_SPA.SparseMatrixCSC)
    for nx in 1:xfmr[:Nx]
        xfmr[:idx_s_fr][nx]   = find_sparse_matrix_idx(y, xfmr[:fr][nx], xfmr[:to][nx])
        xfmr[:idx_s_to][nx]   = find_sparse_matrix_idx(y, xfmr[:to][nx], xfmr[:fr][nx])
        xfmr[:idx_sh_to][nx]  = find_sparse_matrix_idx(y, xfmr[:to][nx], xfmr[:to][nx])
end end

# init #########################################################################
""
function init_admittance_matrix!(branch::NamedTuple, xfmr::NamedTuple, h::Float64)
    # init full matrix idxs i and j, and corresponding values v
    i, j, v = Int[], Int[], ComplexF64[]

    # fill i, j, v for the relevant edges and units
    for nb in 1:branch[:Nb]
        fr = branch[:fr][nb]
        to = branch[:to][nb]
        add!(i, j, v, fr, fr, admittance_branch_series_ondiag(branch[:r][nb], branch[:x][nb], h))
        add!(i, j, v, to, to, admittance_branch_series_ondiag(branch[:r][nb], branch[:x][nb], h))
        add!(i, j, v, fr, to, admittance_branch_series_offdiag(branch[:r][nb], branch[:x][nb], h))
        add!(i, j, v, to, fr, admittance_branch_series_offdiag(branch[:r][nb], branch[:x][nb], h))
        add!(i, j, v, fr, fr, admittance_branch_shunt_ondiag(branch[:b_fr][nb], branch[:g_fr][nb], h))
        add!(i, j, v, to, to, admittance_branch_shunt_ondiag(branch[:b_to][nb], branch[:g_to][nb], h))
    end

    for nx in 1:xfmr[:Nx]
        fr = xfmr[:fr][nx]
        to = xfmr[:to][nx]
        add!(i, j, v, fr, fr, admittance_xfmr_series_ondiag(xfmr[:r][nx], xfmr[:x][nx], h))
        add!(i, j, v, to, to, admittance_xfmr_series_ondiag(xfmr[:r][nx], xfmr[:x][nx], h))
        add!(i, j, v, fr, to, admittance_xfmr_series_offdiag(xfmr[:r][nx], xfmr[:x][nx], h))
        add!(i, j, v, to, fr, admittance_xfmr_series_offdiag(xfmr[:r][nx], xfmr[:x][nx], h))
        add!(i, j, v, to, to, admittance_xfmr_shunt_ondiag(xfmr[:b][nx], xfmr[:g][nx], h))
    end
    
    # create sparse matrix
    y = _SPA.sparse(i, j, v)

    # find idx for the relevant edges and units 
    fill_branch_idx!(branch, y)
    fill_xfmr_idx!(xfmr, y)

    return y
end

# update #######################################################################
""
function update_admittance_matrix!(y::_SPA.SparseMatrixCSC, branch::NamedTuple, xfmr::NamedTuple, nh::Float64)
    # reset the values of the sparse matrix
    y.nzval .= zeros(Complex, length(y.nzval))

    # fill the values of the sparse matrix related to the branches
    for nb in 1:branch[:Nb]
        y[branch[:idx_s_fr][nb]]    += admittance_branch_series_ondiag(branch[:r][nb], branch[:x][nb], nh)
        y[branch[:idx_s_to][nb]]    += admittance_branch_series_ondiag(branch[:r][nb], branch[:x][nb], nh)
        y[branch[:idx_s_fr][nb]]    += admittance_branch_series_offdiag(branch[:r][nb], branch[:x][nb], nh)
        y[branch[:idx_s_to][nb]]    += admittance_branch_series_offdiag(branch[:r][nb], branch[:x][nb], nh)
        y[branch[:idx_sh_fr][nb]]   += admittance_branch_shunt_ondiag(branch[:b_fr][nb], branch[:g_fr][nb], nh)
        y[branch[:idx_sh_to][nb]]   += admittance_branch_shunt_ondiag(branch[:b_to][nb], branch[:g_to][nb], nh)
    end

    # fill the values of the sparse matrix related to the transformers
    for nx in 1:xfmr[:Nx]
        y[xfmr[:idx_s_fr][nx]]     += admittance_xfmr_series_ondiag(xfmr[:r][nx], xfmr[:x][nx], nh)
        y[xfmr[:idx_s_to][nx]]     += admittance_xfmr_series_ondiag(xfmr[:r][nx], xfmr[:x][nx], nh)
        y[xfmr[:idx_s_fr][nx]]     += admittance_xfmr_series_offdiag(xfmr[:r][nx], xfmr[:x][nx], nh)
        y[xfmr[:idx_s_to][nx]]     += admittance_xfmr_series_offdiag(xfmr[:r][nx], xfmr[:x][nx], nh)
        y[xfmr[:idx_sh_to][nx]]    += admittance_xfmr_shunt_ondiag(xfmr[:b][nx], xfmr[:g][nx], nh)
    end

    return y
end

# main #########################################################################
""
function calculate_pos_seq_harmonic_impedance(fdata::Dict{String,Any}, 
                                              h::Vector{Float64}, 
                                              idn::Vector{Int})
    # init harmonic impedance for the relevant nodes
    z       = Dict(ni => Complex[] for ni in idn)

    # init the necessary named tuples for the relevant edges and units
    branch  = init_branch_data(fdata["branch"])
    xfmr    = init_xfmr_data(fdata["xfmr"])

    # init the admittance matrix and update the named tuples
    y       = init_admittance_matrix!(branch, xfmr, h[1])

    # calculate the harmonic impedance for the relevant nodes
    for ni in idn
        i       = zeros(ComplexF64, length(idn))
        i[ni]   = 1.0 + 0.0im

        push!(z[ni], (y \ i)[ni])
    end

    # enumerate over the harmonics
    for nh in h[2:end]
        # update the admittance matrix
        y = update_admittance_matrix!(y, branch, xfmr, nh)
        
        # calculate the harmonic impedance for the relevant nodes
        for ni in idn
            i       = zeros(ComplexF64, length(idn))
            i[ni]   = 1.0 + 0.0im
            push!(z[ni], (y \ i)[ni])
        end 
    end

    return z
end