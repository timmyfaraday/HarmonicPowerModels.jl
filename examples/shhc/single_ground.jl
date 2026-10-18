################################################################################
# Helper: convert a replicated multinetwork to a SINGLE harmonic ground while     #
# keeping the fundamental power flow unchanged ("keep bus 1 as fundamental slack").#
#                                                                               #
#   • bus 1 stays the only reference (bus_type 3) → fundamental = standard PF,     #
#     harmonic ground = bus 1 only.                                              #
#   • the other generator buses become PV (type 2): clean at the fundamental,     #
#     and at HARMONICS they reach ground only through their subtransient          #
#     impedance X'' — added as branches g→1 in the harmonic networks ONLY         #
#     (scaled r''·√h, x''·h), so the fundamental network is never modified.        #
################################################################################

# generator subtransient impedance (IEEE-30 data-file header), fundamental values
const GEN_ZPP = Dict(1=>(0.0025,0.05), 2=>(0.0025,0.05), 5=>(0.025,0.5),
                     8=>(0.025,0.5), 11=>(0.025,0.5), 13=>(0.025,0.5))

function to_single_ground!(hd; gen_buses = collect(keys(GEN_ZPP)), ground = 1,
                           gen_z = GEN_ZPP)
    gset = Set(gen_buses)
    # 1) single reference: ground stays type-3, other gen buses → PV (type 2)
    for (_, ntw) in hd["nw"], (_, b) in ntw["bus"]
        (b["bus_i"] in gset && b["bus_i"] != ground) && (b["bus_type"] = 2)
    end
    # 2) add X'' branches g→ground to HARMONIC networks only (fundamental untouched)
    for (nw, ntw) in hd["nw"]
        nw == "1" && continue
        h = parse(Int, nw)
        tmpl = first(values(ntw["branch"]))          # template for required fields
        nextidx = maximum(b["index"] for b in values(ntw["branch"])) + 1
        for g in gen_buses
            g == ground && continue
            r, x = gen_z[g]
            br = deepcopy(tmpl)
            br["f_bus"] = g; br["t_bus"] = ground
            br["br_r"] = r*sqrt(h); br["br_x"] = x*h
            br["b_fr"] = 0.0; br["b_to"] = 0.0
            br["g_fr"] = 0.0; br["g_to"] = 0.0
            br["tap"]  = 1.0; br["shift"] = 0.0
            br["br_status"] = 1; br["transformer"] = false
            br["index"] = nextidx
            br["f_idx"] = (nextidx, g, ground); br["t_idx"] = (nextidx, ground, g)
            haskey(br,"rate_a") && (br["rate_a"] = 999.0)
            haskey(br,"c_rating") && (br["c_rating"] = 999.0)
            haskey(br,"source_id") && (br["source_id"] = ["branch", nextidx])
            ntw["branch"]["$nextidx"] = br
            nextidx += 1
        end
    end
    return hd
end
