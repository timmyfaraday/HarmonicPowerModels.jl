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

""
init_branch_data(branch::Dict{String,Any}) = 
   (idx_s_fr    = zeros(Int, length(branch)),
    idx_s_to    = zeros(Int, length(branch)),
    idx_sh_fr   = zeros(Int, length(branch)),
    idx_sh_to   = zeros(Int, length(branch)),
    fr          = Int[br["bus_fr"] for br in values(branch)],
    to          = Int[br["bus_to"] for br in values(branch)],
    r           = Float64[br["r"] for br in values(branch)],
    x           = Float64[br["x"] for br in values(branch)],
    b_fr        = Float64[br["b_fr"] for br in values(branch)],
    g_fr        = Float64[br["g_fr"] for br in values(branch)],
    b_to        = Float64[br["b_to"] for br in values(branch)],
    g_to        = Float64[br["g_to"] for br in values(branch)])

""
function init_admittance_matrix_branch!(i::Vector{Int}, j::Vector{Int}, v::Vector{Complex}, h::Float64, branch::NamedTuple)
    for nb in 1:length(branch[:idx])
        # 1) series admittance on-diagonal
        push!(i, branch["fr"][nb])
        push!(j, branch["fr"][nb])
        push!(v, admittance_branch_series_ondiag(branch["r"][nb], branch["x"][nb], h))
        push!(i, branch["to"][nb])
        push!(j, branch["to"][nb])
        push!(v, admittance_branch_series_ondiag(branch["r"][nb], branch["x"][nb], h))
        # 2) series admittance off-diagonal
        push!(i, branch["fr"][nb])
        push!(j, branch["to"][nb])
        push!(v, admittance_branch_series_offdiag(branch["r"][nb], branch["x"][nb], h))
        push!(i, branch["to"][nb])
        push!(j, branch["fr"][nb])
        push!(v, admittance_branch_series_offdiag(branch["r"][nb], branch["x"][nb], h))
        # 3) shunt admittance on-diagonal
        push!(i, branch["fr"][nb])
        push!(j, branch["fr"][nb])
        push!(v, admittance_branch_shunt_ondiag(branch["b_fr"][nb], branch["g_fr"][nb], h))
        push!(i, branch["to"][nb])
        push!(j, branch["to"][nb])
        push!(v, admittance_branch_shunt_ondiag(branch["b_fr"][nb], branch["g_to"][nb], h))
    end
end

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
function fill_branch_idx!(branch::NamedTuple, y::SparseMatrixCSC)
    for nb in 1:length(branch[:idx])
        branch[:idx_s_fr][nb]   = find_sparse_matrix_idx(y, branch[:fr][nb], branch[:to][nb])
        branch[:idx_s_to][nb]   = find_sparse_matrix_idx(y, branch[:to][nb], branch[:fr][nb])
        branch[:idx_sh_fr][nb]  = find_sparse_matrix_idx(y, branch[:fr][nb], branch[:fr][nb])
        branch[:idx_sh_to][nb]  = find_sparse_matrix_idx(y, branch[:to][nb], branch[:to][nb])
end end

""
function find_sparse_index(A,i,j)
    nzr = nzrange(A,j)
    rows = view(rowvals(A),nzr)
    return nzr[searchsortedfirst(rows,i)]
end

""
function init_admittance_matrix!(branch::NamedTuple, nh)
    # init full matrix idxs i and j, and corresponding values v
    i, j, v = Int[], Int[], Complex[]

    # fill i, j, v for the relevant edges and units
    init_admittance_matrix_branch!(i, j, v, branch, nh)    
    # init_admittance_matrix_xfmr!...
    
    # create sparse matrix
    y       = sparse(i, j, v)

    # find idx for the relevant edges and units 
    fill_branch_idx!(branch, y)
    # fill_xfmr_idx!...
end
""
function update_admittance_matrix!(y)
    

end

function calculate_pos_seq_harmonic_impedance(fdata::Dict{String,Any}, 
                                              h::Vector{Float64}, 
                                              idn::Vector{Int})
    # init harmonic impedance for the relevant nodes
    z       = (ni = Complex[] for ni in idn)

    # init the necessary named tuples for the relevant edges and units
    branch  = init_branch_data(fdata["branch"])

    # init the admittance matrix and update the named tuples
    y       = init_admittance_matrix!(branch, h[1])

    # calculate the harmonic impedance for the relevant nodes
    for ni in idn
        i       = zeros(Complex, length(idn))
        i[ni]   = 1.0 + 0.0im
        push!(z[ni], (y / i)[ni])
    end

    # enumerate over the harmonics
    for nh in h[2:end]
        # update the admittance matrix
        y       = update_admittance_matrix!(branch, nh)
        
        # calculate the harmonic impedance for the relevant nodes
        for ni in idn
            i       = zeros(Complex, length(idn))
            i[ni]   = 1.0 + 0.0im
            push!(z[ni], (y / i)[ni])
    end end

    return z
end


function calculate_pos_seq_harmonic_impedance(data::Dict{String,Any}, 
    freq::Vector,
    idn::Vector)
# get necessary parameters
Nn = length(data["bus"])

# init harmonic impedance dictionary
Z = Dict(ni => Complex[] for ni in idn)

# build admittance matrix
Y = build_admittance_matrix(data)

# enumerate of the required frequency
for nf in freq
global h = nf / 50.0

Yh = eval.(Y)

for ni in idn
I = zeros(Complex, Nn)
I[ni] = 1.0 + 0.0im

push!(Z[ni], (Yh \ I)[ni])
end end

# return Z 
return Z
end

# gen
""
function admittance_gen_shunt_ondiag(gen::Dict{String,Any})
    r, x = gen["rsc"], gen["xsc"]

    return :(1 / ($(r) * sqrt(h) + im * $(x) * h))
end

# xfmr
""
function admittance_xfmr_series_ondiag(xfmr::Dict{String,Any})
    r1, r2, x = xfmr["r1"], xfmr["r2"], xfmr["xsc"]

    return :(1 / (sqrt(h) * ($(r1) + $(r2)) + im * h * $(x)))
end
""
function admittance_xfmr_series_offdiag(xfmr::Dict{String,Any})
    r1, r2, x = xfmr["r1"], xfmr["r2"], xfmr["xsc"]

    return :(-1 / (sqrt(h) * ($(r1) + $(r2)) + im * h * $(x)))
end
""
function admittance_xfmr_shunt_ondiag_to(xfmr::Dict{String,Any})
    g   = xfmr["gsh"]

    return :($(g) / sqrt(h))
end

# admittance
@eval function build_admittance_matrix(data::Dict{String,Any})
    # init necessary parameters
    Nn  = length(data["bus"])

    # init admittance matrix 
    Y = [:(0 * h) for ni in 1:Nn, nj in 1:Nn]

    # add the branch admittances
    for (nb, branch) in data["branch"]
        push!(nf,I)
        push!(nt,J)
        
        nf, nt = branch["f_bus"], branch["t_bus"]

        # 1) series admittance on-diagonal
        Y[nf,nf] = :($(Y[nf,nf]) + $(admittance_branch_series_ondiag(branch)))
        Y[nt,nt] = :($(Y[nt,nt]) + $(admittance_branch_series_ondiag(branch)))

        # 2) series admittance off-diagonal
        Y[nf,nt] = :($(Y[nf,nt]) + $(admittance_branch_series_offdiag(branch)))
        Y[nt,nf] = :($(Y[nt,nf]) + $(admittance_branch_series_offdiag(branch)))

        # 3) shunt admittance on-diagonal
        Y[nf,nf] = :($(Y[nf,nf]) + $(admittance_branch_shunt_ondiag_fr(branch)))
        Y[nt,nt] = :($(Y[nt,nt]) + $(admittance_branch_shunt_ondiag_to(branch)))
    end

    # add the xfmr admittances
    if haskey(data, "xfmr")
        for (nx, xfmr) in data["xfmr"]
            nf, nt = xfmr["f_bus"], xfmr["t_bus"]

            # 1) series admittance on-diagonal
            Y[nf,nf] = :($(Y[nf,nf]) + $(admittance_xfmr_series_ondiag(xfmr)))
            Y[nt,nt] = :($(Y[nt,nt]) + $(admittance_xfmr_series_ondiag(xfmr)))

            # 2) series admittance off-diagonal
            Y[nf,nt] = :($(Y[nf,nt]) + $(admittance_xfmr_series_offdiag(xfmr)))
            Y[nt,nf] = :($(Y[nt,nf]) + $(admittance_xfmr_series_offdiag(xfmr)))

            # 3) shunt admittance on-diagonal
            Y[nt,nt] = :($(Y[nt,nt]) + $(admittance_xfmr_shunt_ondiag_to(xfmr)))
        end
    end

    # add the gen admittances
    for (ng, gen) in data["gen"]
        nn  = gen["gen_bus"]
        
        # 1) shunt admittance on-diagonal
        Y[nn,nn] = :($(Y[nn,nn]) + $(admittance_gen_shunt_ondiag(gen)))
    end

    # return admittance matrix
    return Y 
end

# impedance
""
function calculate_pos_seq_harmonic_impedance(data::Dict{String,Any}, 
                                              freq::Vector,
                                              idn::Vector)
    # get necessary parameters
    Nn = length(data["bus"])

    # init harmonic impedance dictionary
    Z = Dict(ni => Complex[] for ni in idn)

    # build admittance matrix
    Y = build_admittance_matrix(data)

    # enumerate of the required frequency
    for nf in freq
        global h = nf / 50.0
        
        Yh = eval.(Y)

        for ni in idn
            I = zeros(Complex, Nn)
            I[ni] = 1.0 + 0.0im

            push!(Z[ni], (Yh \ I)[ni])
    end end

    # return Z 
    return Z
end