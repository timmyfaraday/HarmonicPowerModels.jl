################################################################################
# HarmonicPowerModels.jl                                                       #
# Extension package of PowerModels.jl for Steady-State Power System            #
# Optimization with Power Harmonics.                                           #
# See http://github.com/timmyfaraday/HarmonicPowerModels.jl                    #
################################################################################
# Authors: Tom Van Acker                                                       #
################################################################################
# Changelog:                                                                   #
# v0.2.0 - reviewed TVA                                                        #
# v0.2.1 - reviewed TVA                                                        #
# v0.3.0 - adapted for extended graph representation                           #
################################################################################

""
function ref_add_core!(ref::Dict{Symbol,Any})
    # get the fundamental network ref 
    nf_ref = ref[:it][_PMs.pm_it_sym][:nw][1]

    # enumerate over all harmonic network refs 
    for (nw, nw_ref) in ref[:it][_PMs.pm_it_sym][:nw]
        # conductor ids
        nw_ref[:conductor_ids] = 1:1
        
        # reference buses
        nw_ref[:ref_buses] = Dict{Int,Any}(k => v for (k,v) in nf_ref[:bus] if v["type"] == 3)

        # clean buses
        nw_ref[:clean_buses] = Dict{Int,Any}(k => v for (k,v) in nf_ref[:bus] if v["type"] == 4)

        # branch 
        nw_ref[:arcs_branch_from]   = [(nb,br["bus"]...)            for (nb,br) in nf_ref[:branch]]
        nw_ref[:arcs_branch_to]     = [(nb,reverse(br["bus"])...)   for (nb,br) in nf_ref[:branch]]
        nw_ref[:arcs_branch]        = [nw_ref[:arcs_branch_from]; nw_ref[:arcs_branch_to]]        
        
        bus_arcs_branch = Dict((i, []) for (i,bus) in nw_ref[:bus])
        for (b,i,j) in nw_ref[:arcs_branch]
            push!(bus_arcs_branch[i], (b,i,j))
        end
        nw_ref[:bus_arcs_branch] = bus_arcs_branch

        # xfmr 
        nw_ref[:arcs_xfmr_from]     = [(nx,xf["bus"]...)            for (nx,xf) in nf_ref[:xfmr]]
        nw_ref[:arcs_xfmr_to]       = [(nx,reverse(xf["bus"])...)   for (nx,xf) in nf_ref[:xfmr]]
        nw_ref[:arcs_xfmr]          = [nw_ref[:arcs_xfmr_from]; nw_ref[:arcs_xfmr_to]]

        bus_arcs_xfmr = Dict((i, []) for (i,bus) in nw_ref[:bus])
        for (x,i,j) in nw_ref[:arcs_xfmr]
            push!(bus_arcs_xfmr[i], (x,i,j))
        end
        nw_ref[:bus_arcs_xfmr] = bus_arcs_xfmr

        # gen
        bus_gen = Dict((i, Int[]) for (i,bus) in nw_ref[:bus])
        for (g,gen) in nw_ref[:gen]
            push!(bus_gen[nf_ref[:gen][g]["bus"]], g)
        end
        nw_ref[:bus_gen] = bus_gen

        # hload
        bus_hload = Dict((i, Int[]) for (i,bus) in nw_ref[:bus])
        for (l,hload) in nw_ref[:hload]
            push!(bus_hload[nf_ref[:hload][l]["bus"]], l)
        end
        nw_ref[:bus_hload] = bus_hload

        # hsrc
        bus_hsrc = Dict((i, Int[]) for (i,bus) in nw_ref[:bus])
        for (s,hsrc) in nw_ref[:hsrc]
            push!(bus_hsrc[nf_ref[:hsrc][s]["bus"]], s)
        end
        nw_ref[:bus_hsrc] = bus_hsrc
    end
end