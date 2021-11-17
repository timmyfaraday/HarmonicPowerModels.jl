################################################################################
#  Copyright 2021, Frederik Geth, Tom Van Acker                                #
################################################################################
# HarmonicPowerModels.jl                                                       #
# An extention package of PowerModels.jl for Harmonic (Optimal) Power Flow     #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################

"""
    HarmonicPowerModels.magnetic_flux_density_polar(E::Vector{<:Real}, θ::Vector{<:Real}, ω::Vector{<:Real}, t::Vector{<:Real}, A::Real, Vbase::Real)

Function to determine the time-domain magnetic flux density B(t) for a given
time-vector `t`, primary turns `N` and core surface `A`, based on the frequency-
domain excitation voltage in polar form, given by `E`, `θ`, `ω` and `Vbase`.

```math 
\begin{align}
    B(t) = \sum_{h \in H} \frac{Vbase ⋅ |E_h|}{A ⋅ N ⋅ ω_h} ⋅ \cos(\omega_h ⋅ t + \theta_h)
\end{align}
```

"""
magnetic_flux_density_polar(E::Vector{<:Real}, θ::Vector{<:Real}, ω::Vector{<:Real}, 
                            t::Vector{<:Real}, A::Real, N::Int, Vbase::Real) = 
    sum(Vbase .* E[h] ./ ω[h] ./ A ./ N .* cos.(ω[h] .* t .+ θ[h]) for h in 1:length(E))

"""
    HarmonicPowerModels.magnetic_flux_density_rectangular(Ere::Vector{<:Real}, Eim::Vector{<:Real}, ω::Vector{<:Real}, t::Vector{<:Real}, A::Real, Vbase::Real)

Function to determine the time-domain magnetic flux density B(t) for a given
time-vector `t`, primary turns `N` and core surface `A`, based on the frequency-
domain excitation voltage in rectangular form, given by `Ere`, `Eim`, `ω` and `Vbase`. 

This dispatches to `magnetic_flux_density_polar(hypot.(Ere,Eim), atan.(Eim,Eim), ω, t, A, N, Vbase)`.
"""
magnetic_flux_density_rectangular(Ere::Vector{<:Real}, Eim::Vector{<:Real}, ω::Vector{<:Real}, 
                                  t::Vector{<:Real}, A::Real, N::Int, Vbase::Real) =
    magnetic_flux_density_polar(hypot.(Ere, Eim), atan.(Eim, Ere), ω, t, A, N, Vbase)

""
function sample_excitation_voltage_polar(data, nx, exc)
    nb = data["nw"]["1"]["xfmr"]["$nx"]["fr_bus"]
    θ_min = [0.0 for (nw,nh) in data["harmonics"] if nh in exc["Hᴱ"]]
    θ_max = [2pi for (nw,nh) in data["harmonics"] if nh in exc["Hᴱ"]]
    dθ    = (θ_max - θ_min) / 9
    E_min = [data["nw"]["$nw"]["bus"]["$nb"]["vmin"] for (nw,nh) in data["harmonics"]
                                                     if nh in exc["Hᴱ"]]
    E_max = [data["nw"]["$nw"]["bus"]["$nb"]["vmax"] for (nw,nh) in data["harmonics"]
                                                     if nh in exc["Hᴱ"]]
    dE    = (E_max - E_min) / 9

    S = reduce(vcat,[[E_min[ni]:dE:E_max[ni],θ_min[ni]:dθ:θ_max[ni]] 
                      for (ni,nh) in enumerate(exc["Hᴱ"])])
    R = [1:length(s) for s in S]
    return S, R
end

""
function sample_voltage_rectangular(data, nx, exc)
    nb = data["nw"]["1"]["xfmr"]["$nx"]["fr_bus"]
    Ere_min = [-data["nw"]["$nw"]["bus"]["$nb"]["vmax"] for (nw,nh) in data["harmonics"]
                                                        if nh in exc["Hᴱ"]]
    Ere_max = [ data["nw"]["$nw"]["bus"]["$nb"]["vmax"] for (nw,nh) in data["harmonics"]
                                                        if nh in exc["Hᴱ"]]
    dEre    = (Ere_max - Ere_min) / 9
    Eim_min = [-data["nw"]["$nw"]["bus"]["$nb"]["vmax"] for (nw,nh) in data["harmonics"]
                                                        if nh in exc["Hᴱ"]]
    Eim_max = [ data["nw"]["$nw"]["bus"]["$nb"]["vmax"] for (nw,nh) in data["harmonics"]
                                                        if nh in exc["Hᴱ"]]
    dEim    = (Eim_max - Eim_min) / 9

    S = reduce(vcat,[[Ere_min[ni]:dEre:Ere_max[ni],Eim_min[ni]:dEim:Eim_max[ni]] 
                      for (ni,nh) in enumerate(E_harmonics)])
    R = [1:length(s) for s in S]
    return S, R
end

"""
    HarmonicPowerModels.sample_xfmr_excitation(data::Dict{String, <:Any}, xfmr_exc::Dict{Int, Dict{String, <:Any})

This function creates anonymous functions which wrap a spline model of the 
exitation current, either in `:polar` or `:rectangular` coordinates. 
As inputs it takes excitation voltage, either in `:polar` or 
`:rectangular` coordinates, of transformer x ∈ 𝓧:
    `$E^{re}_{h,x}$, $E^{im}_{h,x}$, ∀ h ∈ 𝓗ᵉ,`
        or
    `$E_{h,x}$, $θ_{h,x}$, ∀ h ∈ 𝓗ᵉ,`
and outputs either the `:polar` or `:rectangular` coordinates of the exictation 
current of tranformer x ∈ 𝓧:
    `$I^{exc,re}_{h,x}$, $I^{exc,im}_{h,x}$, ∀ h ∈ 𝓗ⁱ,`
        or
    `$I^{exc}_{h,x}$, $φ^{exc}_{h,x}$, ∀ h ∈ 𝓗ⁱ,`
where 𝓗ᵉ and 𝓗ⁱ denote the set of excitation voltage and current harmonics, 
respectively.
"""
function sample_xfmr_magnetizing_current(data::Dict{String, <:Any}, xfmr_exc::Dict{Int, Dict{String, <:Any}})
    # set of all harmonics and corresponding nw ids
    idx = sortperm(collect(values(data["harmonics"])))
    NW  = collect(keys(data["harmonics"]))[idx]                                 # sorted set of corresponding networks [String]
    H   = collect(values(data["harmonics"]))[idx]                               # sorted set of harmonics [Int]

    for (nx, exc) in enumerate(xfmr_exc)
        # assert all necessary keys are in exc
        @assert isempty(setdiff(["Hᴱ", "Hᴵ" "Fᴱ", "Fᴵ", "l", "A", "N", "BH"], keys(exc)))

        # define the decomposition structure, see https://github.com/JuliaDynamics/SignalDecomposition.jl 
        fq = _SDC.Sinusoidal(_HPM.freq .* H)

        # get the base voltage
        nb    = data["nw"]["1"]["xfmr"][nx]["f_bus"]
        Vbase = data["nw"]["1"]["bus"][nb]["base_kv"] * 1e3
        Abase = data["nw"]["1"]["baseMVA"] * 1e6 / Vbase

        # determine the necessary time parameters
        ωᴱ = (2.0 * pi * _HPM.freq) .* exc["Hᴱ"]                                # angular frequency [rad/Hz] for all relevant excitation voltage harmonics
        dt = (1 / (100 * _HPM.freq * maximum(exc["Hᴵ"])))                       # time-step of the time-domain excitation voltage, length of full wave of highest considered current harmonic divided by 100
        t  = 0.0:dt:(5.0 / _HPM.freq)                                           # time-range of the time-domain excitation voltage

        # sample the excitation voltage
        if exc["Fᴱ"] == :polar
            S, R = sample_excitation_voltage_polar(data, nx, exc)
        elseif exc["Fᴱ"] == :rectangular
            S, R = sample_excitation_voltage_rectangular(data, nx, exc)
        else
            error("Fᴱ ∉ [:polar, :rectangular] for xfmr $nx")
        end
        
        # initialize dictionaries for excitation current, where:
        # I_formulation == :polar => a = magnitude, b = phase angle
        # I_formulation == :rectangular => a = real, b = imaginary
        Ia = Dict(nh => zeros(R...) for nh in exc["Hᴵ"])
        Ib = Dict(nh => zeros(R...) for nh in exc["Hᴵ"])

        # sample the magnetizing current
        @showprogress for nr in Iterators.product(R...)
            # get a excitation voltage sample
            sample = [S[ni][ns] for (ni,ns) in enumerate(nr)]

            # determine the time-domain magnetic flux density B(t) [T] based on 
            # the excitation voltage sample
            if exc["Fᴱ"] == :polar
                E, θ = sample[1:2:end], sample[2:2:end]
                B = magnetic_flux_density_polar(E, θ, ωᴱ, t, exc["A"], exc["N"], Vbase)
            else exc["Fᴱ"] == :rectangular
                Ere, Eim = sample[1:2:end], sample[2:2:end]
                B = magnetic_flux_density_rectangular(Ere, Eim, ωᴱ, t, exc["A"], exc["N"], Vbase)
            end

            # determine the time-domain magnetizing current iᵐ(t) [pu] based on 
            # the magnetic field intensity H(t) = BH(B(t)) [A-turn/m]
            Im = exc["l"] / exc["N"] .* exc["BH"].(B) ./ Abase

            # decompose the time-domain magnetizing current iᵐ(t) [pu] in its
            # frequency-domain components for all harmonics
            _SDC.decompose(t, Im, fq)

            # translate the frequency components to the required magnetizing 
            # current formulation
            if exc["Fᴵ"] == :polar
                # NOTE -- angle convention is reversed -> introduce minus-sign 
                # for the phase angle
                I, φ = fq.A[2:end], -fq.φ[2:end]
                for (ni,nh) in enumerate(H)
                    Ia[nh][nr...], Ib[nh][nr...] = I[ni], φ[ni]
                end
            elseif exc["Fᴵ"] == :rectangular
                # NOTE -- angle convention is reversed -> introduce minus-sign 
                # for the phase angle
                Ire, Iim = fq.A[2:end] .* sin.(-fq.φ[2:end]), fq.A[2:end] .* cos.(-fq.φ[2:end])
                for (ni,nh) in enumerate(H)
                    Ia[nh][nr...], Ib[nh][nr...] = Ire[ni], Iim[ni]
                end 
            else
                error("Fᴵ ∉ [:polar, :rectangular] for xfmr $nx")
            end 
        end

        # fill the xfmr data structure, enumerating over all harmonics 
        for (nw,nh) in data["harmonics"]
            # shortcut for the xfmr data
            bus  = data["nw"][nw]["bus"]["$nb"]
            xfmr = data["nw"][nw]["xfmr"]["$nx"]

            # set general data
            xfmr["Fᴱ"]  = exc["Fᴱ"]
            xfmr["Fᴵ"]  = exc["Fᴵ"]
            xfmr["NWᴱ"] = NW[[nh in exc["Hᵉ"] for nh in H]]
            xfmr["NWᴵ"] = NW[[nh in exc["Hᴵ"] for nh in H]]

            # interpolate and set magnetizing current data
            if nh in exc["Hᴵ"]
                method = _INT.BSpline(_INT.Cubic(_INT.Line(_INT.OnGrid())))
                xfmr["INT_A"] = _INT.extrapolate(_INT.scale(_INT.interpolate(Ia[nh], method), S...), _INT.Line())
                xfmr["INT_B"] = _INT.extrapolate(_INT.scale(_INT.interpolate(Ia[nh], method), S...), _INT.Line())
                xfmr["Im_A"]  = (x...) -> xfmr["INT_A"](x...)
                xfmr["Im_B"]  = (x...) -> xfmr["INT_B"](x...)
                xfmr["dIm_A"] = (x...) -> _INT.gradient(xfmr["INT_A"], x...)
                xfmr["dIm_B"] = (x...) -> _INT.gradient(xfmr["INT_A"], x...)
            end

            # set the excitation voltage limits
            xfmr["eat_min"], xfmr["eat_max"] = 0.0, 2π
            xfmr["emt_min"], xfmr["emt_max"] = bus["vmin"], bus["vmax"]
            xfmr["ert_min"], xfmr["ert_max"] = -bus["vmax"], bus["vmax"]
            xfmr["eit_min"], xfmr["eit_max"] = -bus["vmax"], bus["vmax"]
        end
    end
end