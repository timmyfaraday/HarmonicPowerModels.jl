################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Authors: Tom Van Acker                                                       #
################################################################################
# Changelog:                                                                   #
# v0.2.2 - init                                                                #
################################################################################


# The idea is to build a function 'build_admittance_matrix() that builds 
# a expression 'Y()', the latter of should allow for efficient evaluation for a 
# given frequency, here expressed relative to the fundamental, i.e., as harmonic
# h. Note that h does not have to be an integer

# This leads us to the idea of metaprogramming, i.e., where evaluation of the 
# frequency only happens at runtime.

# Currently, we consider branches, transformers and generators

# Note that for any considered aspect here, a numerical value is needed, given 
# that a zero-value, i.e., short-circuit, in an admittance-matrix is Inf.

# The unit of all admittances, and therefore resulting impedances is pu, and 
# should be transformed using the relevant base for the respective voltage level

# A conversation with J. Tant is advised to improve the tractability of the 
# current implementation. As defining h as a global will not be most efficient.

# branch, pi-model #############################################################
## series component
### on-diagonal, i.e., (i,i)
#### zs     = sqrt(h) * r + im * h * x 
#### ys     = 1 / (sqrt(h) * r + im * h * x)
#### yᵢᵢ    = 1 / (sqrt(h) * r + im * h * x)
function admittance_branch_series_ondiag(branch::Dict{String,Any})
    r, x = branch["br_r"], branch["br_x"]

    return :(1 / (sqrt(h) * $(r) + im * h * $(x)))
end

### off-diagonal, i.e., (i,j)
#### zs     = sqrt(h) * r + im * h * x 
#### ys     = 1 / (sqrt(h) * r + im * h * x)
#### yᵢⱼ    = -1 / (sqrt(h) * r + im * h * x)
function admittance_branch_series_offdiag(branch::Dict{String,Any})
    r, x = branch["br_r"], branch["br_x"]

    return :(-1 / (sqrt(h) * $(r) + im * h * $(x)))
end

## shunt component, on-diagonal, fr-side
#### ysh    = im * b / h
#### yᵢᵢ    = im * b / h
function admittance_branch_shunt_ondiag_fr(branch::Dict{String,Any})
    b   = branch["b_fr"]

    return :(im * $(b) * h)
end

## shunt component, on-diagonal, to-side
#### ysh    = im * b / h
#### yᵢᵢ    = im * b / h
function admittance_branch_shunt_ondiag_to(branch::Dict{String,Any})
    b   = branch["b_to"]

    return :(im * $(b) * h)
end

# transformer, shifted t-model, i.e., shunt at to node #########################
## series component
### on-diagonal, i.e., (i,i)
#### zs     = sqrt(h) * (r₁ + r₂) + im * h * xˢᶜ 
#### ys     = 1 / (sqrt(h) * (r₁ + r₂) + im * h * xˢᶜ)
#### yᵢᵢ    = 1 / (sqrt(h) * (r₁ + r₂) + im * h * xˢᶜ)
function admittance_xfmr_series_ondiag(xfmr::Dict{String,Any})
    r1, r2, x = xfmr["r1"], xfmr["r2"], xfmr["xsc"]

    return :(1 / (sqrt(h) * ($(r1) + $(r2)) + im * h * $(x)))
end

### off-diagonal, i.e., (i,j)
#### zs     = sqrt(h) * (r₁ + r₂) + im * h * x 
#### ys     = 1 / (sqrt(h) * (r₁ + r₂) + im * h * x)
#### yᵢⱼ    = -1 / (sqrt(h) * (r₁ + r₂) + im * h * x)
function admittance_xfmr_series_offdiag(xfmr::Dict{String,Any})
    r1, r2, x = xfmr["r1"], xfmr["r2"], xfmr["xsc"]

    return :(-1 / (sqrt(h) * ($(r1) + $(r2)) + im * h * $(x)))
end

## shunt component, on-diagonal, to-side
#### ysh    = g / sqrt(h)
#### yᵢᵢ    = g / sqrt(h)
function admittance_xfmr_shunt_ondiag_to(xfmr::Dict{String,Any})
    g   = xfmr["gsh"]

    return :($(g) / sqrt(h))
end

# generator, shunt at gen node #################################################
#### ysh    = 1 / (rˢᶜ * sqrt(h) + im * h * xˢᶜ)
#### yᵢᵢ    = 1 / (rˢᶜ * sqrt(h) + im * h * xˢᶜ)
function admittance_gen_shunt_ondiag(gen::Dict{String,Any})
    r, x = gen["rsc"], gen["xsc"]

    return :(1 / ($(r) * sqrt(h) + im * $(x) * h))
end

# motor, shunt at mot node #####################################################
#### ysh    = 1 / (rˢᶜ * sqrt(h) + im * h * xˢᶜ)
#### yᵢᵢ    = 1 / (rˢᶜ * sqrt(h) + im * h * xˢᶜ)
function admittance_mot_shunt_ondiag(mot::Dict{String,Any})
    r, x = mot["rsc"], mot["xsc"]

    return :(1 / ($(r) * sqrt(h) + im * $(x) * h))
end

# capacitance, shunt at mot node #####################################################
#### ysh    = 1 / (rˢᶜ * sqrt(h) + im * h * xˢᶜ)
#### yᵢᵢ    = 1 / (rˢᶜ * sqrt(h) + im * h * xˢᶜ)
function admittance_cap_shunt_ondiag(cap::Dict{String,Any})
    xl, xc = cap["xl"], cap["xc"]

    return :(1 / (im * $(xl) * h - im * $(xc) / h))
end

# inf bus, shunt at inf node ###################################################
function admittance_inf_shunt_ondiag(inf::Dict{String,Any})
    H, z, a = inf["h"], inf["z"], inf["a"]

    r = _INT.linear_interpolation(H, z .* cos.(a))
    x = _INT.linear_interpolation(H, z .* sin.(a))

    return :(1 / ($r(h) + im * $x(h)))
end

# build admittance matrix ######################################################
@eval function build_admittance_matrix(data::Dict{String,Any})
    # init necessary parameters
    Nn  = length(data["bus"])

    # init admittance matrix 
    Y = [:(0 * h) for ni in 1:Nn, nj in 1:Nn]

    # add the branch admittances
    for (nb, branch) in data["branch"]
        nf_, nt_ = branch["f_bus"], branch["t_bus"]

        nf = data["bus"]["$nf_"]["hb_idx"]
        nt = data["bus"]["$nt_"]["hb_idx"]

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
            # nf, nt = xfmr["f_bus"], xfmr["t_bus"]

            nf_, nt_ = xfmr["f_bus"], xfmr["t_bus"]

            nf = data["bus"]["$nf_"]["hb_idx"]
            nt = data["bus"]["$nt_"]["hb_idx"]

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
        nn_  = gen["gen_bus"]
        nn = data["bus"]["$nn_"]["hb_idx"]

        # 1) shunt admittance on-diagonal
        Y[nn,nn] = :($(Y[nn,nn]) + $(admittance_gen_shunt_ondiag(gen)))
    end

    # add motor admittances
    for (nm, mot) in data["mot"]
        nn_ = mot["mot_bus"]
        nn  = data["bus"]["$nn_"]["hb_idx"]

        # 1) shunt admittance on-diagonal
        Y[nn,nn] = :($(Y[nn,nn]) + $(admittance_mot_shunt_ondiag(mot)))
    end

    # add capacitance admittances
    for (nc, cap) in data["cap"]
        nn_ = cap["cap_bus"]
        nn  = data["bus"]["$nn_"]["hb_idx"]

        # 1) shunt admittance on-diagonal
        Y[nn,nn] = :($(Y[nn,nn]) + $(admittance_cap_shunt_ondiag(cap)))
    end

    # add the inf admittances
    for (ni, inf) in data["inf"]
        nn_  = inf["inf_bus"]
        nn   = data["bus"]["$nn_"]["hb_idx"]

        # 1) shunt admittance on-diagonal
        Y[nn,nn] = :($(Y[nn,nn]) + $(admittance_inf_shunt_ondiag(inf)))
    end

    # return admittance matrix
    return Y 
end

global h = 1

# calculate harmonic impedance
function calculate_pos_seq_harmonic_impedance(data::Dict{String,Any}, 
                                              freq::Vector,
                                              idn::Vector)
    # init necessary parameters
    Nn  = length(data["bus"])
    
    # init harmonic impedance dictionary
    Z = Dict(ni => Complex[] for ni in idn)

    # build admittance matrix

    println("Building Y matrix")

    @time Y = build_admittance_matrix(data)

    # enumerate of the required frequency
    for nf in freq
        global h = nf / 50.0

        println("Evaluate Yh for ", nf, " Hz")
        @time Yh = eval.(Y)

        Yhc = zeros(ComplexF64, Nn, Nn)

        for i in 1:Nn
            for j in 1:Nn
                Yhc[i, j] = float(Yh[i,j])
            end
        end

        Yhs = _SPA.sparse(Yhc)

        for ni in idn
            I = zeros(ComplexF64, Nn)
            I[ni] = 1.0 + 0.0im

            println("Calculate Z for bus ", ni)
            @time push!(Z[ni], (Yhs \ I)[ni])
    end end
 
    return Z
end