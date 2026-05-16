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
#### ysh    = rˢᶜ * sqrt(h) + im * h * xˢᶜ
#### yᵢᵢ    = rˢᶜ * sqrt(h) + im * h * xˢᶜ
function admittance_gen_shunt_ondiag(gen::Dict{String,Any})
    r, x = gen["rsc"], gen["xsc"]

    return :(1 / ($(r) * sqrt(h) + im * $(x) * h))
end

# build admittance matrix ######################################################
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

    # add the gen admittances
    for (ng, gen) in data["gen"]
        nn  = gen["gen_bus"]
        
        # 1) shunt admittance on-diagonal
        Y[nn,nn] = :($(Y[nn,nn]) + $(admittance_gen_shunt_ondiag(gen)))
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

# small-case sanity tests ######################################################
## single branch (r = 0.01, x=0.1, b / 2 = 0.01)
f = 750.0

r = 0.01
x = 0.1
b = 0.01

### method
data = Dict{String,Any}("bus"       => Dict{String,Any}("$n" => [] for n in 1:2),
                        "branch"    => Dict{String,Any}("1" => Dict{String,Any}("f_bus" => 1,
                                                                                "t_bus" => 2,
                                                                                "br_r"  => r,
                                                                                "br_x"  => x,
                                                                                "b_fr"  => b,
                                                                                "b_to"  => b)),
                        "gen"       => Dict{String,Any}(),
                        "xfmr"      => Dict{String,Any}())
Zc = calculate_pos_seq_harmonic_impedance(data, [f], [1,2])
### hand-calc
R = r
L = x / 2 / pi / 50
C = b / 2 / pi / 50

Zl = R * sqrt(f / 50) + im * 2 * pi * f * L + 1 / (im * 2 * pi * f * C)
Zr = 1 / (im * 2 * pi * f * C)

Zm = 1 / (1 / Zl + 1 / Zr)

### tests 
@assert Zc[1][1] ≈ Zm
@assert Zc[1][1] ≈ Zc[2][1]

## single generator (Ssc=2.1GVA, XR=20) #######################################
f = 425.0

Ssc     = 2.1e9
Sbase   = 100e6
XRr     = 20.0

z = Sbase / Ssc
r = z / sqrt(1 + XRr^2)
x = sqrt(z^2 - r^2)

@assert x == z / sqrt(1 + (1/XRr)^2)

### method
data = Dict{String,Any}("bus"       => Dict{String,Any}("$n" => [] for n in 1),
                        "branch"    => Dict{String,Any}(),
                        "gen"       => Dict{String,Any}("1" => Dict{String,Any}("gen_bus" => 1,
                                                                                "rsc" => r,
                                                                                "xsc" => x)),
                        "xfmr"      => Dict{String,Any}())
Zc = calculate_pos_seq_harmonic_impedance(data, [f], [1])

### hand-calc
R = r
L = x / 2 / pi / 50

Zm = R * sqrt(f / 50) + im * 2 * pi * f * L

### tests 
@assert Zc[1][1] ≈ Zm

## single transformer ##########################################################
f = 695.0

x = 0.208
g = 0.00001
r = 2 * 0.004988662

### method
data = Dict{String,Any}("bus"       => Dict{String,Any}("$n" => [] for n in 1:2),
                        "branch"    => Dict{String,Any}(),
                        "gen"       => Dict{String,Any}(),
                        "xfmr"      => Dict{String,Any}("1" => Dict{String,Any}("f_bus" => 1,
                                                                                "t_bus" => 2,
                                                                                "r1" => r / 2,
                                                                                "r2" => r / 2,
                                                                                "xsc" => x,
                                                                                "gsh" => g)))
Zc = calculate_pos_seq_harmonic_impedance(data, [f], [1,2])

### hand-calc, from node 1
R = r + 1 / g
L = x / 2 / pi / 50

Zm = R * sqrt(f / 50.0) + im * 2 * pi * f * L 

### tests 
@assert Zc[1][1] ≈ Zm

## 30bus System ################################################################
using HarmonicPowerModels
using PowerModels
using Plots

const PMs = PowerModels
const HPM = HarmonicPowerModels

path = joinpath(HPM.BASE_DIR,"test/data/matpower/30bus_network_hhc.m")
data = PMs.parse_file(path)

Ssc     = 2.1e9
Sbase   = 100e6
XRr     = 20.0

z = Sbase / Ssc
r = z / sqrt(1 + XRr^2)
x = sqrt(z^2 - r^2)
for (ng, gen) in data["gen"]
    gen["rsc"] = r
    gen["xsc"] = x
end

Zc = calculate_pos_seq_harmonic_impedance(data, collect(50.0:10.0:2500.0), collect(1:30))

# plot 
plot(50.0:10.0:2500.0, abs.(Zc[9]), yaxis=:log)