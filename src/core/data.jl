"""
    HarmonicPowerModels.replicate(data::Dict{String,<:Any};
                                  harmonics::Array{Int}=Int[],
                                  xfmr_exc::Bool=false)

Creates a multinetwork data-file for all considered harmonics based on the data
specified for the fundamental harmonic.

The considered harmonics are provided by the user through the named 
argument harmonics, and collected based on the harmonic bus data.

Additionally, a Bool `xfmr_exc` may be passed to include the excitation current
of the transformers.

NOTES:
i) branch
-   The susceptance b of any branch is assumed to be purely capacitive, i.e., 
    b = ω * C.
-   The resistance r of any branch is scaled using the square of the harmonic.
ii) transformer
-   The transformer configuration defaults to Yy0.
-   The resistance r of any transformer is scaled using the square of the harmonic.
"""
function _HPM.replicate(data::Dict{String, Any}; 
                        harmonics::Array{Int}=Int[],
                        xfmr_exc::Dict{String, Any}=Dict{String, Any}())
    # extend the user-provided harmonics based on the data
    collect_harmonics!(data, harmonics, xfmr_exc)

    # create a multinetwork data structure
    Nh = length(harmonics)
    data = _PMs.replicate(data, Nh)
    
    # add the harmonics to the overall data
    data["harmonics"] = Dict("$nw" => harmonics[nw] for nw in 1:Nh)

    # add the xfmr excitation
    sample_xfmr_excitation(data, xfmr_exc)

    # re-evaluate the data for each harmonic 
    for nw in keys(data["nw"])
        nh = data["harmonics"][nw]
        mh = 2 * pi * nh

        # re-evaluate load power
        for load in values(data["nw"][nw]["load"])
            nb = load["source_id"][2]
            bus = data["nw"][nw]["bus"]["$nb"]
            multiplier = haskey(bus, "nh_$nh") ? bus["nh_$nh"] : 0.0 ;
            haskey(load, "pd") ? load["pd"] *= multiplier : ~ ;
            haskey(load, "qd") ? load["qd"] *= multiplier : ~ ;
        end

        # re-evaluate the bus data 
        for bus in values(data["nw"][nw]["bus"])
            #use fundamental as limit for rms
            bus["vminrms"] = bus["vmin"]
            bus["vmaxrms"] = bus["vmax"]
            
            # true harmonics don't have minimum voltage
            if nw !=1 
                bus["vmin"] = 0 
            end
        end

        # re-evaluate the branch data 
        for branch in values(data["nw"][nw]["branch"])
            haskey(branch, "br_r") ? branch["br_r"] *= sqrt(nh) : ~ ;
            haskey(branch, "br_x") ? branch["br_x"] *= nh : ~ ;
            haskey(branch, "b_fr") ? branch["b_fr"] *= nh : ~ ;
            haskey(branch, "b_to") ? branch["b_to"] *= nh : ~ ;
        end

        # re-evaluate the transformer data
        if haskey(data["nw"][nw], "xfmr")
            for xfmr in values(data["nw"][nw]["xfmr"])
                haskey(xfmr, "xsc") ? xfmr["xsc"] *= mh : ~ ;
                haskey(xfmr, "r1")  ? xfmr["r1"] *= sqrt(nh) : ~ ;
                haskey(xfmr, "r2")  ? xfmr["r2"] *= sqrt(nh) : ~ ;

                haskey(xfmr, "xe1") ? xfmr["xe1"] *= mh : ~ ;
                haskey(xfmr, "xe2") ? xfmr["xe2"] *= mh : ~ ;
                haskey(xfmr, "re1") ? xfmr["re1"] *= sqrt(nh) : ~ ;
                haskey(xfmr, "re2") ? xfmr["re2"] *= sqrt(nh) : ~ ;

                xfmr["cnf1"] = haskey(xfmr, "vg") ? uppercase(xfmr["vg"][1]) : "Y" ;
                xfmr["cnf2"] = haskey(xfmr, "vg") ? uppercase(xfmr["vg"][2]) : "Y" ;

                shift = haskey(xfmr, "vg") ? parse(Int, xfmr["vg"][3]) : 0 ;
                if is_pos_sequence(nh)
                    xfmr["tr"] = cosd(-30.0 * shift)
                    xfmr["ti"] = sind(-30.0 * shift)
                elseif is_neg_sequence(nh)
                    xfmr["tr"] = cosd(30.0 * shift)
                    xfmr["ti"] = sind(30.0 * shift)
                elseif is_zero_sequence(nh)
                    xfmr["tr"] = 1.0
                    xfmr["ti"] = 0.0
                end
            end
        end
    end

    return data
end

function collect_harmonics!(data::Dict{String, Any}, harmonics::Array{Int}, xfmr_exc::Dict{String, Any})
    bus_harmonics = [parse(Int,ni[4:end])   for ni in keys(data["bus"]["1"]) 
                                            if  ni[1:2] == "nh"]
    push!(harmonics, bus_harmonics...)
    if haskey(xfmr_exc, "voltage_harmonics") 
        push!(harmonics, xfmr_exc["voltage_harmonics"]...)
    end
    if haskey(xfmr_exc, "current_harmonics")
        push!(harmonics, xfmr_exc["current_harmonics"]...)
    end  
    unique!(sort!(harmonics))
end

is_pos_sequence(nh::Int) = nh % 3 == 1
is_neg_sequence(nh::Int) = nh % 3 == 2
is_zero_sequence(nh::Int) = nh % 3 == 0


### TRANSFORMER EXCITATION ###
""
excitation_flux_polar(V, θ, w, t) = 
    sum(V[n] ./ w[n] .* sin.(w[n] .* t .+ θ[n]) for n in 1:length(V))
excitation_flux_rectangular(Vre, Vim, w, t) =
    excitation_flux_polar(sqrt.(Vre.^2 .+ Vim.^2), atan.(Vim./(Vre.+1e-8)), w, t)
excitation_current_sigmoid(inom, ψmax, ψ) = 
    -inom .* log.(2.0 ./ (ψ ./ ψmax .+ 1.0) .- 1.0)

""
function sample_voltage_polar(voltage_harmonics, dv, vmin, vmax, dθ, θmin, θmax)
    S = reduce(vcat,[[vmin[ni]:dv[ni]:vmax[ni],θmin[ni]:dθ[ni]:θmax[ni]] 
                      for (ni,nh) in enumerate(voltage_harmonics)])
    R = [1:length(s) for s in S]
    return S, R
end
function sample_voltage_rectangular(voltage_harmonics, dv, vmin, vmax)
    S = reduce(vcat,[[vmin[ni]:dv[ni]:vmax[ni],vmin[ni]:dv[ni]:vmax[ni]] 
                      for (ni,nh) in enumerate(voltage_harmonics)])
    R = [1:length(s) for s in S]
    return S, R
end

""
function sample_xfmr_excitation(data::Dict{String, <:Any}, xfmr_exc::Dict{String, <:Any})
    if isempty(xfmr_exc)
        ## HARMONICS
        voltage_harmonics = []
        current_harmonics = []
    else
        ## HARMONICS
        voltage_harmonics = xfmr_exc["voltage_harmonics"]
        current_harmonics = xfmr_exc["current_harmonics"]

        ## TIME
        N = xfmr_exc["N"]

        ## EXCITATION CURRENT 
        current_type = xfmr_exc["current_type"]

        ## EXCITATION FLUX 
        excitation_type = xfmr_exc["excitation_type"] 
        inom = xfmr_exc["inom"]
        ψmax = xfmr_exc["ψmax"]

        ## EXCITATION VOLTAGE 
        voltage_type = xfmr_exc["voltage_type"]
        dv   = xfmr_exc["dv"]
        vmin = xfmr_exc["vmin"]
        vmax = xfmr_exc["vmax"]
        dθ   = xfmr_exc["dθ"]
        θmin = xfmr_exc["θmin"]
        θmax = xfmr_exc["θmax"]

        ## TIME
        w = (2.0 * pi * freq) .* voltage_harmonics
        t = 0.0:(1 / (freq * N * maximum(current_harmonics))):(5.0 / freq)

        ## DECOMPOSITION 
        fq = _SDC.Sinusoidal(freq .* current_harmonics)

        ## SAMPLE VOLTAGE
        if voltage_type == :polar
            S, R = sample_voltage_polar(voltage_harmonics, dv, vmin, vmax, dθ, θmin, θmax)
        elseif voltage_type == :rectangular
            S, R = sample_voltage_rectangular(voltage_harmonics, dv, vmin, vmax)
        end

        ## INITIALIZE EXCITATION CURRENT DICTIONARIES
        Ia = Dict(nh => zeros(R...) for nh in current_harmonics)
        Ib = Dict(nh => zeros(R...) for nh in current_harmonics)

        ## SAMPLING LOOP   
        @showprogress for nr in Iterators.product(R...)
            sample = [S[ni][ns] for (ni,ns) in enumerate(nr)]

            if voltage_type == :polar
                V, θ = sample[1:2:end], sample[2:2:end]
                ψexc = excitation_flux_polar(V, θ, w, t)
            else voltage_type == :rectangular 
                Vre, Vim = sample[1:2:end], sample[2:2:end]
                ψexc = excitation_flux_rectangular(Vre, Vim, w, t) 
            end

            if excitation_type == :sigmoid
                I_exc = excitation_current_sigmoid(inom, ψmax, ψexc)
            end

            _SDC.decompose(t, I_exc, fq)

            if current_type == :polar 
                I, φ = fq.A[2:end], fq.φ[2:end]
                for (ni,nh) in enumerate(current_harmonics)
                    Ia[nh][nr...], Ib[nh][nr...] = I[ni], φ[ni]
                end
            elseif current_type == :rectangular
                Ire, Iim = fq.A[2:end] .* sin.(fq.φ[2:end]), fq.A[2:end] .* cos.(fq.φ[2:end])
                for (ni,nh) in enumerate(current_harmonics)
                    Ia[nh][nr...], Ib[nh][nr...] = Ire[ni], Iim[ni]
                end
            end
        end
    end
    reverse_harmonics = Dict(value => key for (key, value) in data["harmonics"])
    current_harmonics_ntws = 
        [reverse_harmonics[nc] for nc in current_harmonics]
    voltage_harmonics_ntws = 
        [reverse_harmonics[nv] for nv in voltage_harmonics]

    method = _INT.BSpline(_INT.Cubic(_INT.Line(_INT.OnGrid())))
    for nw in keys(data["nw"]) 
        ni = parse(Int, nw)
        nh = data["harmonics"][nw]
        if haskey(data["nw"][nw], "xfmr")
            for xfmr in values(data["nw"][nw]["xfmr"])
                xfmr["current_harmonics_ntws"] = current_harmonics_ntws
                xfmr["voltage_harmonics_ntws"] = voltage_harmonics_ntws
                if nh in current_harmonics
                    xfmr["EXC_A"]  = _INT.scale(_INT.interpolate(Ia[nh], method), S...)
                    xfmr["EXC_B"]  = _INT.scale(_INT.interpolate(Ib[nh], method), S...)
                    xfmr["INT_A"]  = (x...) -> xfmr["EXC_A"](x...)
                    xfmr["INT_B"]  = (x...) -> xfmr["EXC_B"](x...)
                    xfmr["GRAD_A"] = (x...) -> _INT.gradient(xfmr["EXC_A"], x...)
                    xfmr["GRAD_B"] = (x...) -> _INT.gradient(xfmr["EXC_B"], x...)
                end
                if nh in voltage_harmonics
                    xfmr["ert_min"], xfmr["ert_max"] = vmin[ni], vmax[ni] 
                    xfmr["eit_min"], xfmr["eit_max"] = vmin[ni], vmax[ni]
                else
                    xfmr["ert_min"], xfmr["ert_max"] = 0.0, 1.1 
                    xfmr["eit_min"], xfmr["eit_max"] = 0.0, 1.1
                end
            end 
        end
    end
end

