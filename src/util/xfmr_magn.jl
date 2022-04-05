################################################################################
#  Copyright 2021, Frederik Geth, Tom Van Acker                                #
################################################################################
# HarmonicPowerModels.jl                                                       #
# An extention package of PowerModels.jl for Harmonic (Optimal) Power Flow     #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################

""
vmin(data, nw, nb) = data["nw"]["$nw"]["bus"]["$nb"]["vmin"]
vmax(data, nw, nb) = data["nw"]["$nw"]["bus"]["$nb"]["vmax"]

"""
    HarmonicPowerModels.sample_excitation_voltage(nx::Int, data{String, <:Any}, 
                                                  xfmr_magn::Dict{String, <:Any})
    
"""
function sample_excitation_voltage(nx::Int, data::Dict{String, <:Any}, xfmr_magn::Dict{String, <:Any})
    # determine the from bus of the transformer nx
    nb = data["nw"]["1"]["xfmr"]["$nx"]["f_bus"]

    # given polar coordinates, determine the min and max of |E| and θ
    if xfmr_magn["Fᴱ"] == :polar
        Ea_min = [vmin(data, nw, nb) for (nw,nh) in data["harmonics"]
                                     if nh in xfmr_magn["Hᴱ"]]
        Ea_max = [vmax(data, nw, nb) for (nw,nh) in data["harmonics"]
                                     if nh in xfmr_magn["Hᴱ"]]
        Eb_min = [0.0 for (nw,nh) in data["harmonics"] if nh in xfmr_magn["Hᴱ"]]
        Eb_max = [2pi for (nw,nh) in data["harmonics"] if nh in xfmr_magn["Hᴱ"]]
    end

    # given rectangular coordinates, determine the min and max of Ere and Eim
    if xfmr_magn["Fᴱ"] == :rectangular 
        Ea_min = [-vmax(data, nw, nb) for (nw,nh) in data["harmonics"]
                                      if nh in xfmr_magn["Hᴱ"]]
        Ea_max = [vmax(data, nw, nb) for (nw,nh) in data["harmonics"]
                                     if nh in xfmr_magn["Hᴱ"]]
        Eb_min = [-vmax(data, nw, nb) for (nw,nh) in data["harmonics"]
                                      if nh in xfmr_magn["Hᴱ"]]
        Eb_max = [vmax(data, nw, nb) for (nw,nh) in data["harmonics"]
                                     if nh in xfmr_magn["Hᴱ"]]
    end

    S = reduce(vcat,[[range(Ea_min[ni], Ea_max[ni], length=10),
                      range(Eb_min[ni], Eb_max[ni], length=10)] 
                      for (ni,nh) in enumerate(xfmr_magn["Hᴱ"])])
    R = [1:length(s) for s in S]
    
    return S, R
end

"""
    HarmonicPowerModels.magnetizing_current(Ea::Vector{<:Real}, Eb::Vector{<:Real},
                                            data::Dict{String, <:Any}, xfmr::Dict{String, <:Any})

"""
function magnetizing_current(Ea::Vector{<:Real}, Eb::Vector{<:Real}, 
                             magn_data::Dict{String, <:Any}, xfmr::Dict{String, <:Any})
    
    # assert all necessary keys are in data and xfmr
    # NOTE: This may be deleted at a later stage once code is streamlined
    @assert isempty(setdiff(["Abase", "fq", "Fᴱ", "Fᴵ", "Hᴱ", "Vbase", "t"], keys(magn_data)))
    @assert isempty(setdiff(["A", "BH", "l", "N"], keys(xfmr)))

    # given a rectangular excitation voltage [pu, pu], translate them to polar coordinates [pu, (rad)]
    # Ea --> real parts, Eb --> imag parts
    if magn_data["Fᴱ"] == :rectangular
        Ea = hypot.(Ea, Eb) 
        Eb = atan.(Eb, Ea)
    end

    # determine the magnetizing current in the time domain [pu]
    # B(t) = √2 * Vbase / (A ⋅ N) * ∑ |Eₕ| / wₕ ⋅ cos(wₕ ⋅ t + ∠Eₕ)
    # E, θ [pu, (rad)] → B(t) [T] → H(t) [A-t/m] → im(t) [pu]
    ω   = 2 .* pi .* _HPM.freq .* magn_data["Hᴱ"]
    B   = sum((sqrt(2) * magn_data["Vbase"] * Ea[ni]) ./ (xfmr["A"] * xfmr["N"] * ω[ni]) .* cos.(ω[ni] .* magn_data["t"] .+ Eb[ni]) for (ni,nh) in enumerate(magn_data["Hᴱ"]))
    Im  = xfmr["l"] ./ (magn_data["Abase"] * xfmr["N"]) .* xfmr["BH"].(B)

    # decompose and return the magnetizing current in the frequency domain [pu, (rad)]
    _SDC.decompose(magn_data["t"], Im, magn_data["fq"])

    # return the magnetizing current in the frequency domain in the required coordinates
    # NOTE -- angle convention is reversed -> introduce minus-sign for the phase angle
    # NOTE -- sqrt(2) convert from amplitude to rms magnitude
    magn_data["Fᴵ"] == :polar && return magn_data["fq"].A[2:end] ./ sqrt(2), 
                                        -magn_data["fq"].φ[2:end]
    magn_data["Fᴵ"] == :rectangular && return   magn_data["fq"].A[2:end] ./ sqrt(2) .* sin.(-magn_data["fq"].φ[2:end]), 
                                                magn_data["fq"].A[2:end] ./ sqrt(2) .* cos.(-magn_data["fq"].φ[2:end])
end

"""""
    HarmonicPowerModels.sample_magnetizing_current(data::Dict{String, <:Any}, 
                                                   xfmr_exc::Dict{String, <:Any})


"""
function sample_magnetizing_current(data::Dict{String,<:Any}, xfmr_magn::Dict{String,<:Any})
    # assert all necessary keys are in xfrm_magn
    @assert isempty(setdiff(["Hᴱ", "Hᴵ", "Fᴱ", "Fᴵ", "xfmr"], keys(xfmr_magn)))
    
    # set of all harmonics and corresponding nw ids
    idx = sortperm(collect(values(data["harmonics"])))                          # sorted set of corresponding networks [String]
    H   = collect(values(data["harmonics"]))[idx]
    NW  = collect(keys(data["harmonics"]))[idx]                                 # sorted set of harmonics [Int]

    # set a data dictionary for decomposition of the magnetizing current, where:
    # - fq: decompostion structure, see https://github.com/JuliaDynamics/SignalDecomposition.jl 
    # - t: timeseries with time-step equal to length of full wave of highest considered current harmonic divided by 100
    magn_data = Dict("Abase" => 0.0,
                     "Fᴱ"    => xfmr_magn["Fᴱ"],
                     "Fᴵ"    => xfmr_magn["Fᴵ"],
                     "fq"    => _SDC.Sinusoidal(_HPM.freq .* H),
                     "Hᴱ"    => xfmr_magn["Hᴱ"],
                     "Vbase" => 0.0,
                     "t"     => 0.0 : (1 / (100 * _HPM.freq * maximum(H))) : (5.0 / _HPM.freq))

    for (nx, xfmr) in xfmr_magn["xfmr"]
        # set the base voltage and current for the xfmr
        nb                  = data["nw"]["1"]["xfmr"]["$nx"]["f_bus"]
        magn_data["Vbase"]  = data["nw"]["1"]["bus"]["$nb"]["base_kv"] * 1e3
        magn_data["Abase"]  = data["nw"]["1"]["baseMVA"] * 1e6 / magn_data["Vbase"]

        # sample the excitation voltage
        S, R = sample_excitation_voltage(nx, data, xfmr_magn)
        
        # initialize dictionaries for magnitizing current, where:
        # I_formulation == :polar => a = magnitude, b = phase angle
        # I_formulation == :rectangular => a = real, b = imaginary
        Ia = Dict(nh => zeros(R...) for nh in xfmr_magn["Hᴵ"])
        Ib = Dict(nh => zeros(R...) for nh in xfmr_magn["Hᴵ"])

        # sample the magnetizing current
        @showprogress for nr in Iterators.product(R...)
            # get a excitation voltage sample
            sample = [S[ni][ns] for (ni,ns) in enumerate(nr)]

            # determine the magnitizing current in the required coordinates
            Ix, Iy = magnetizing_current(sample[1:2:end], sample[2:2:end], magn_data, xfmr)
            
            # fill the magnitizing current dictionaries Ia and Ib
            for (ni,nh) in enumerate(H)
                Ia[nh][nr...], Ib[nh][nr...] = Ix[ni], Iy[ni]
            end 
        end

        # fill the xfmr data structure, enumerating over all harmonics 
        for (nw,nh) in data["harmonics"]
            # shortcut for the xfmr data
            bus  = data["nw"][nw]["bus"]["$nb"]
            xfrm = data["nw"][nw]["xfmr"]["$nx"] # note xfrm ≠ xfmr

            # set general data
            xfrm["Fᴱ"]  = xfmr_magn["Fᴱ"]
            xfrm["Fᴵ"]  = xfmr_magn["Fᴵ"]
            xfrm["NWᴱ"] = NW[[nh in xfmr_magn["Hᴱ"] for nh in H]]
            xfrm["NWᴵ"] = NW[[nh in xfmr_magn["Hᴵ"] for nh in H]]

            # interpolate and set magnetizing current data
            if nh in H
                method = _INT.BSpline(_INT.Cubic(_INT.Line(_INT.OnGrid())))
                xfrm["INT_A"] = _INT.extrapolate(_INT.scale(_INT.interpolate(Ia[nh], method), S...), _INT.Line())
                xfrm["INT_B"] = _INT.extrapolate(_INT.scale(_INT.interpolate(Ia[nh], method), S...), _INT.Line())
                xfrm["Im_A"]  = (x...) -> xfrm["INT_A"](x...)
                xfrm["Im_B"]  = (x...) -> xfrm["INT_B"](x...)
                xfrm["dIm_A"] = (x...) -> _INT.gradient(xfrm["INT_A"], x...)
                xfrm["dIm_B"] = (x...) -> _INT.gradient(xfrm["INT_A"], x...)
            end

            # set the excitation voltage limits
            xfrm["eat_min"], xfrm["eat_max"] = 0.0, 2π
            xfrm["emt_min"], xfrm["emt_max"] = bus["vmin"], bus["vmax"]
            xfrm["ert_min"], xfrm["ert_max"] = -bus["vmax"], bus["vmax"]
            xfrm["eit_min"], xfrm["eit_max"] = -bus["vmax"], bus["vmax"]
        end
    end
end