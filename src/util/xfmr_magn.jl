################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Authors: Tom Van Acker, Frederik Geth                                        #
################################################################################
# Changelog:                                                                   #
# v0.2.0 - reviewed TVA                                                        #
################################################################################

"""
    sample_magnetizing_current
"""
function sample_magnetizing_current(hdata::Dict{String,<:Any}, xfmr_magn::Dict{String,<:Any})
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
        Abase   = hdata["nw"]["1"]["baseMVA"] * 10e6 / xfmr["Vbase"]            # base current [A]

        # sample the excitation voltage
        IDH, Emax, pcs = xfmr_magn["IDH"], xfmr_magn["Emax"], xfmr_magn["pcs"]
        S = [range(-IDH[ni] * Emax, IDH[ni] * Emax, length=pcs[ni]) 
                for ni in repeat(1:NHᴱ, inner=2)]                               # samples of the real and imaginary excitation voltage, consecutively, for each excitation voltage harmonic numbers
        R = [1:pcs[ni] for ni in repeat(1:NHᴱ, inner=2)]                        # range of the samples of the real and imaginary excitation voltage, consecutively, for each excitation voltage harmonic numbers

        # initialize dictionaries for magnitizing current, where:
        # magnitizing current output
        Ire = Dict(nw => zeros(repeat(pcs, inner=2)...) for nw in keys(hdata["nw"]))
        Iim = Dict(nw => zeros(repeat(pcs, inner=2)...) for nw in keys(hdata["nw"]))

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
            for (ni,nh) in enumerate(keys(hdata["nw"]))
                Ire[nh][nr...] = Imgn[ni] .* cos(Iphs[ni])                      # determine the real part [pu] of the magnitizing current, consecutively, for each magnitizing current harmonic number
                Iim[nh][nr...] = Imgn[ni] .* sin(Iphs[ni])                      # determine the imaginary part [pu] of the magnitizing current, consecutively, for each magnitizing current harmonic number
            end
        end

        # fill the xfmr data structure, enumerating over all harmonics 
        for (nw,ntw) in hdata["nw"]
            # shortcut for the xfmr data
            bus  = ntw["bus"]["$(ntw["xfmr"]["$nx"]["f_bus"])"]
            xfrm = ntw["xfmr"]["$nx"]                                           # note xfrm ≠ xfmr

            # set general data
            xfrm["Hᴱ"] = xfmr_magn["Hᴱ"]
            xfrm["Hᴵ"] = xfmr_magn["Hᴵ"]

            # interpolate and set magnetizing current data
            method = _INT.BSpline(_INT.Cubic(_INT.Line(_INT.OnGrid())))
            xfrm["INT_A"] = _INT.extrapolate(_INT.scale(_INT.interpolate(Ire[nw], method), S...), _INT.Line())
            xfrm["INT_B"] = _INT.extrapolate(_INT.scale(_INT.interpolate(Iim[nw], method), S...), _INT.Line())
            xfrm["Im_A"]  = (x...) -> xfrm["INT_A"](x...)
            xfrm["Im_B"]  = (x...) -> xfrm["INT_B"](x...)

            # set the excitation voltage limits
            xfrm["eax_min"], xfrm["eax_max"] = 0.0, 2π
            xfrm["emx_min"], xfrm["emx_max"] = bus["vmin"], bus["vmax"]
            xfrm["erx_min"], xfrm["erx_max"] = - xfmr_magn["Emax"], xfmr_magn["Emax"]
            xfrm["eix_min"], xfrm["eix_max"] = - xfmr_magn["Emax"], xfmr_magn["Emax"]
        end
    end
end