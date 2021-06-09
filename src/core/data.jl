"""
    HarmonicPowerModels.replicate(data::Dict{String,<:Any})

Creates a multinetwork data-file for all considered harmonics based on the data
specified for the fundamental harmonic.

The considered harmonics are provided by the user through the named 
argument harmonics, and collected based on the harmonic bus data.

NOTES:
i) branch
-   The susceptance b of any branch is assumed to be purely capacitive, i.e., 
    b = ω * C.
-   The resistance r of any branch is scaled using the square of the harmonic.
ii) transformer
-   The transformer configuration defaults to Yy0.
-   The resistance r of any transformer is scaled using the square of the harmonic.
"""
function replicate(data::Dict{String, <:Any}; harmonics::Array{Int}=Int[])
    # extend the user-provided harmonics based on the data
    collect_harmonics!(data, harmonics)

    # create a multinetwork data structure
    Nh = length(harmonics)
    data = _PMs.replicate(data, Nh)
    
    # add the harmonics to the overall data
    data["harmonics"] = Dict("$nw" => harmonics[nw] for nw in 1:Nh)

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

        # re-evaluate the bus data -- NEEDS DISCUSSION
        # for bus in values(data["nw"][nw]["bus"])

        # end

        # re-evaluate the branch data 
        for branch in values(data["nw"][nw]["branch"])
            haskey(branch, "x") ? branch["x"] *= mh : ~ ;
            haskey(branch, "b") ? branch["b"] *= mh : ~ ;
            haskey(branch, "r") ? branch["r"] *= sqrt(nh) : ~ ;
        end

        # re-evaluate the transformer data
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

    return data
end

function collect_harmonics!(data::Dict{String, <:Any}, harmonics::Array{Int})
    bus_harmonics = [parse(Int,ni[4:end])   for ni in keys(data["bus"]["1"]) 
                                            if  ni[1:2] == "nh"]
    push!(harmonics, bus_harmonics...)
    unique!(sort!(harmonics))
end


is_pos_sequence(nh::Int) = nh % 3 == 1
is_neg_sequence(nh::Int) = nh % 3 == 2
is_zero_sequence(nh::Int) = nh % 3 == 0

