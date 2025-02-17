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

# branch 
""
function admittance_branch_series_ondiag(branch::Dict{String,Any})
    r, x = branch["br_r"], branch["br_x"]

    return :(1 / (sqrt(h) * $(r) + im * h * $(x)))
end
""
function admittance_branch_series_offdiag(branch::Dict{String,Any})
    r, x = branch["br_r"], branch["br_x"]

    return :(-1 / (sqrt(h) * $(r) + im * h * $(x)))
end
""
function admittance_branch_shunt_ondiag_fr(branch::Dict{String,Any})
    b   = branch["b_fr"]

    return :(im * $(b) * h)
end
""
function admittance_branch_shunt_ondiag_to(branch::Dict{String,Any})
    b   = branch["b_to"]

    return :(im * $(b) * h)
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