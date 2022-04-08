################################################################################
#  Copyright 2021, Frederik Geth, Tom Van Acker                                #
################################################################################
# HarmonicPowerModels.jl                                                       #
# An extention package of PowerModels.jl for Harmonic (Optimal) Power Flow     #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################

# ""
# vmin(data, nw, nb) = data["nw"]["$nw"]["bus"]["$nb"]["vmin"]
# vmax(data, nw, nb) = data["nw"]["$nw"]["bus"]["$nb"]["vmax"]

# """
#     HarmonicPowerModels.sample_excitation_voltage(nx::Int, data{String, <:Any}, 
#                                                   xfmr_magn::Dict{String, <:Any})
    
# """
# function sample_excitation_voltage(nx::Int, data::Dict{String, <:Any}, xfmr_magn::Dict{String, <:Any})
#     # determine the from bus of the transformer nx
#     nb = data["nw"]["1"]["xfmr"]["$nx"]["f_bus"]

#     # given polar coordinates, determine the min and max of |E| and θ
#     if xfmr_magn["Fᴱ"] == :polar
#         Ea_min = [vmin(data, nw, nb) for (nw,nh) in data["harmonics"]
#                                      if nh in xfmr_magn["Hᴱ"]]
#         Ea_max = [vmax(data, nw, nb) for (nw,nh) in data["harmonics"]
#                                      if nh in xfmr_magn["Hᴱ"]]
#         Eb_min = [0.0 for (nw,nh) in data["harmonics"] if nh in xfmr_magn["Hᴱ"]]
#         Eb_max = [2pi for (nw,nh) in data["harmonics"] if nh in xfmr_magn["Hᴱ"]]
#     end

#     # given rectangular coordinates, determine the min and max of Ere and Eim
#     if xfmr_magn["Fᴱ"] == :rectangular 
#         Ea_min = [-vmax(data, nw, nb) for (nw,nh) in data["harmonics"]
#                                       if nh in xfmr_magn["Hᴱ"]]
#         Ea_max = [vmax(data, nw, nb) for (nw,nh) in data["harmonics"]
#                                      if nh in xfmr_magn["Hᴱ"]]
#         Eb_min = [-vmax(data, nw, nb) for (nw,nh) in data["harmonics"]
#                                       if nh in xfmr_magn["Hᴱ"]]
#         Eb_max = [vmax(data, nw, nb) for (nw,nh) in data["harmonics"]
#                                      if nh in xfmr_magn["Hᴱ"]]
#     end

#     S = reduce(vcat,[[range(Ea_min[ni], Ea_max[ni], length=10),
#                       range(Eb_min[ni], Eb_max[ni], length=10)] 
#                       for (ni,nh) in enumerate(xfmr_magn["Hᴱ"])])
#     R = [1:length(s) for s in S]
    
#     return S, R
# end

# """
#     HarmonicPowerModels.magnetizing_current(Ea::Vector{<:Real}, Eb::Vector{<:Real},
#                                             data::Dict{String, <:Any}, xfmr::Dict{String, <:Any})

# """
# function magnetizing_current(Ea::Vector{<:Real}, Eb::Vector{<:Real}, 
#                              magn_data::Dict{String, <:Any}, xfmr::Dict{String, <:Any})
    
#     # assert all necessary keys are in data and xfmr
#     # NOTE: This may be deleted at a later stage once code is streamlined
#     @assert isempty(setdiff(["Abase", "fq", "Fᴱ", "Fᴵ", "Hᴱ", "Vbase", "t"], keys(magn_data)))
#     @assert isempty(setdiff(["A", "BH", "l", "N"], keys(xfmr)))

#     # given a rectangular excitation voltage [pu, pu], translate them to polar coordinates [pu, (rad)]
#     # Ea --> real parts, Eb --> imag parts
#     if magn_data["Fᴱ"] == :rectangular
#         Ea = hypot.(Ea, Eb) 
#         Eb = atan.(Eb, Ea)
#     end

#     # determine the magnetizing current in the time domain [pu]
#     # B(t) = √2 * Vbase / (A ⋅ N) * ∑ |Eₕ| / wₕ ⋅ cos(wₕ ⋅ t + ∠Eₕ)
#     # E, θ [pu, (rad)] → B(t) [T] → H(t) [A-t/m] → im(t) [pu]
#     ω   = 2 .* pi .* _HPM.freq .* magn_data["Hᴱ"]
#     B   = sum((sqrt(2) * magn_data["Vbase"] * Ea[ni]) ./ (xfmr["A"] * xfmr["N"] * ω[ni]) .* cos.(ω[ni] .* magn_data["t"] .+ Eb[ni]) for (ni,nh) in enumerate(magn_data["Hᴱ"]))
#     Im  = xfmr["l"] ./ (magn_data["Abase"] * xfmr["N"]) .* xfmr["BH"].(B)

#     # decompose and return the magnetizing current in the frequency domain [pu, (rad)]
#     _SDC.decompose(magn_data["t"], Im, magn_data["fq"])

#     # return the magnetizing current in the frequency domain in the required coordinates
#     # NOTE -- angle convention is reversed -> introduce minus-sign for the phase angle
#     # NOTE -- sqrt(2) convert from amplitude to rms magnitude
#     magn_data["Fᴵ"] == :polar && return magn_data["fq"].A[2:end] ./ sqrt(2), 
#                                         -magn_data["fq"].φ[2:end]
#     magn_data["Fᴵ"] == :rectangular && return   magn_data["fq"].A[2:end] ./ sqrt(2) .* sin.(-magn_data["fq"].φ[2:end]), 
#                                                 magn_data["fq"].A[2:end] ./ sqrt(2) .* cos.(-magn_data["fq"].φ[2:end])
# end

"""""
    HarmonicPowerModels.sample_magnetizing_current(data::Dict{String, <:Any}, 
                                                   xfmr_exc::Dict{String, <:Any})


"""
function sample_magnetizing_current(data::Dict{String,<:Any}, xfmr_magn::Dict{String,<:Any})
    # set of all harmonics and corresponding nw ids
    idx = sortperm(collect(values(data["harmonics"])))                          # sorted set of corresponding networks [String]
    HW  = collect(values(data["harmonics"]))[idx]
    NW  = collect(keys(data["harmonics"]))[idx]                                 # sorted set of harmonics [Int]

    # derived input
    dt      = (1 / (100 * _HPM.freq * maximum(xfmr_magn["Hᴵ"])))
    tmax    = (5.0 / _HPM.freq)
    t       = 0.0 : dt : tmax                                                   # time range [s]
    ωᴱ      = 2.0 * pi * _HPM.freq .* xfmr_magn["Hᴱ"]                           # angular frequencies of the excitation voltage harmonic numbers [(rad)/s]
    NHᴱ     = length(xfmr_magn["Hᴱ"])                                           # number of excitation voltage harmonic numbers [-]
    NHᴵ     = length(xfmr_magn["Hᴵ"])                                           # number of magnitizing current harmonic numbers [-]

    # decomposition input 
    fq      = _SDC.Sinusoidal(_HPM.freq .* xfmr_magn["Hᴵ"])

    for (nx, xfmr) in xfmr_magn["xfmr"]
        Abase   = data["nw"]["1"]["baseMVA"] / xfmr["Vbase"]                    # base current [A]

        # sample the excitation voltage
        IDH, Emax, pcs = xfmr_magn["IDH"], xfmr_magn["Emax"], xfmr_magn["pcs"]
        S = [range(-IDH[ni] * Emax, IDH[ni] * Emax, length=pcs[ni]) 
                for ni in repeat(1:NHᴱ, inner=2)]                               # samples of the real and imaginary excitation voltage, consecutively, for each excitation voltage harmonic numbers
        R = [1:pcs[ni] for ni in repeat(1:NHᴱ, inner=2)]                        # range of the samples of the real and imaginary excitation voltage, consecutively, for each excitation voltage harmonic numbers

        # initialize dictionaries for magnitizing current, where:
        # magnitizing current output
        Ire = Dict(nh => zeros(repeat(pcs, inner=2)...) for nh in HW)
        Iim = Dict(nh => zeros(repeat(pcs, inner=2)...) for nh in HW)

        # sample the magnetizing current
        @showprogress for nr in Iterators.product(R...)
            # get a excitation voltage sample
            sample = [S[ni][ns] for (ni,ns) in enumerate(nr)]

            Ere = xfmr["Vbase"] / sqrt(3) .* sample[1:2:end]                    # select the real part [V] of the excitation voltage sample, consecutively, for each excitation voltage harmonic number
            Eim = xfmr["Vbase"] / sqrt(3) .* sample[2:2:end]                    # select the imaginary part [V] of the excitation voltage sample, consecutively, for each excitation voltage harmonic number
            
            Emgn = sqrt(2) .* hypot.(Ere, Eim)                                  # create the magnitude [V] of the excitation voltage sample, consecutively, for each excitation voltage harmonic number
            Ephs = atan.(Eim, Ere)                                              # create the phase of the excitation voltage sample, consecutively, for each excitation voltage harmonic number

            # E(t) = ∑ₕ |Eₕ| ⋅ sin(ωₕ ⋅ t + θₕ)                                          
            # NOTE: THE SINUS CONVENTION IS USED
            E   = sum(Emgn[ni] .* sin.(ωᴱ[ni] .* t .+ Ephs[ni]) for ni in 1:NHᴱ)
            # B(t) = 1 / (A ⋅ N) ∫ -E(t) dt
            #      = 1 / (A ⋅ N) ∑ₕ |Eₕ| / ωₕ ⋅ cos(ωₕ ⋅ t + θₕ)     
            B   = 1 / (xfmr["A"] * xfmr["N"]) .* sum(Emgn[ni] ./ ωᴱ[ni] .* 
                                cos.(ωᴱ[ni] .* t .+ Ephs[ni]) for ni in 1:NHᴱ)  # determine the magnetic flux density [T] based on the excitation voltage
            # H(t) = fᴮᴴ(B(t))
            H   = xfmr["BH"].(B)                                                # determine the magnitic field intensity [A/m] based on the BH curve
            # I(t) = l / N ⋅ H(t)
            Im  = xfmr["l"] / xfmr["N"] .* H                                    # determine the magnetizing current [A]

            # decompose the magnetizing current and fill the matrices
            _SDC.decompose(t, Im, fq)
            Imgn = fq.A[2:end] ./ (sqrt(2) * Abase)                             # determine the magnitude [pu] of the magnitizing current, consecutively, for each magnitizing current harmonic number
            # NOTE: THE DECOMPOSE RETURNS COSINUS, WE USE THE SINUS CONVENTION,
            # HENCE THE ADDITION OF PI / 2 
            Iphs = fq.φ[2:end] .+ pi / 2                                        # determine the phase [(rad)] of the magnitizing current, consecutively, for each magnitizing current harmonic number
            for (ni,nh) in enumerate(HW)
                Ire[nh][nr...] = Imgn[ni] .* cos(Iphs[ni])                      # determine the real part [pu] of the magnitizing current, consecutively, for each magnitizing current harmonic number
                Iim[nh][nr...] = Imgn[ni] .* sin(Iphs[ni])                      # determine the imaginary part [pu] of the magnitizing current, consecutively, for each magnitizing current harmonic number
            end
        end

        # fill the xfmr data structure, enumerating over all harmonics 
        for (nw,nh) in data["harmonics"]
            # shortcut for the xfmr data
            nb   = data["nw"]["1"]["xfmr"]["$nx"]["f_bus"]
            bus  = data["nw"][nw]["bus"]["$nb"]
            xfrm = data["nw"][nw]["xfmr"]["$nx"] # note xfrm ≠ xfmr

            # set general data
            xfrm["Fᴱ"]  = xfmr_magn["Fᴱ"]
            xfrm["Fᴵ"]  = xfmr_magn["Fᴵ"]
            xfrm["NWᴱ"] = NW[[nh in xfmr_magn["Hᴱ"] for nh in HW]]
            xfrm["NWᴵ"] = NW[[nh in xfmr_magn["Hᴵ"] for nh in HW]]

            # interpolate and set magnetizing current data
            if nh in HW
                method = _INT.BSpline(_INT.Cubic(_INT.Line(_INT.OnGrid())))
                xfrm["INT_A"] = _INT.extrapolate(_INT.scale(_INT.interpolate(Ire[nh], method), S...), _INT.Line())
                xfrm["INT_B"] = _INT.extrapolate(_INT.scale(_INT.interpolate(Iim[nh], method), S...), _INT.Line())
                xfrm["Im_A"]  = (x...) -> xfrm["INT_A"](x...)
                xfrm["Im_B"]  = (x...) -> xfrm["INT_B"](x...)
                xfrm["dIm_A"] = (x...) -> _INT.gradient(xfrm["INT_A"], x...)
                xfrm["dIm_B"] = (x...) -> _INT.gradient(xfrm["INT_A"], x...)
            end

            # set the excitation voltage limits
            xfrm["eat_min"], xfrm["eat_max"] = 0.0, 2π
            xfrm["emt_min"], xfrm["emt_max"] = bus["vmin"], bus["vmax"]
            xfrm["ert_min"], xfrm["ert_max"] = - xfmr_magn["Emax"], xfmr_magn["Emax"]
            xfrm["eit_min"], xfrm["eit_max"] = - xfmr_magn["Emax"], xfmr_magn["Emax"]
        end
    end
end