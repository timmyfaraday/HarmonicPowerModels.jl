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

    # for (t, xfmr) in data["xfmr"]
    #     xfmr["voltage_harmonics"] = []
    #     xfmr["current_harmonics"] = []
    #     xfmr["voltage_harmonics_ntws"] = []
    #     xfmr["current_harmonics_ntws"] = []
    # end
    @show harmonics
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
            load["multiplier"] = multiplier
        end

        # re-evaluate gen 
        for gen in values(data["nw"][nw]["gen"])
            if haskey(gen, "isfilter") && gen["isfilter"] == 1
                # do nothing, is handled by constraints constraint_active_filter
                gen["pmin"] = -abs(gen["pmax"])
                gen["qmin"] = -abs(gen["qmax"])
            else #is true generator
                if nw !="1" #cost of harmonics set to 0 
                    gen["cost"] *= 0 
                    #harmonics can be injected/absorbed to match load 
                    gen["pmin"] = -abs(gen["pmax"])
                    gen["qmin"] = -abs(gen["qmax"])
                end
            end

        end

        # re-evaluate the bus data 
        for bus in values(data["nw"][nw]["bus"])
            #use fundamental as limit for rms
            bus["vminrms"] = bus["vmin"]
            bus["vmaxrms"] = bus["vmax"]
            
            # true harmonics don't have minimum voltage
            if nw !="1" 
                bus["vmin"] = 0 
            end
            rm = Dict(3=> 0.05, 5=>0.06, 7=> 0.05, 9=> 0.015, 11=>0.035)
            # inject relative magnitude limits for harmonics
            bus["rm"] = rm
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
                haskey(xfmr, "xsc") ? xfmr["xsc"] *= nh : ~ ;
                haskey(xfmr, "r1")  ? xfmr["r1"] *= sqrt(nh) : ~ ;
                haskey(xfmr, "r2")  ? xfmr["r2"] *= sqrt(nh) : ~ ;

                haskey(xfmr, "xe1") ? xfmr["xe1"] *= nh : ~ ;
                haskey(xfmr, "xe2") ? xfmr["xe2"] *= nh : ~ ;
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
