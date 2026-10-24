################################################################################
# HarmonicPowerModels.jl — shared setup for the SHHC paper scripts              #
#                                                                               #
# The paper scripts (paper_equivalence.jl, paper_capacitive.jl,                  #
# paper_types.jl, paper_bounds.jl) each produce the data of one subsection of   #
# the numerical illustration, as CSV files in output/paper/. All use the IEEE    #
# 30-bus system with a purely inductive base network (no line charging), the    #
# local angle reference and N = 21 harmonic units unless stated otherwise.     #
# Author: Tom Van Acker                                                          #
################################################################################

include(joinpath(@__DIR__, "recomb_ieee30.jl"))
include(joinpath(@__DIR__, "converter_types.jl"))

const NS        = 21
const PAPER_OUT = joinpath(@__DIR__, "output", "paper"); mkpath(PAPER_OUT)
const PAPER_SCENARIOS = ["inductive", "inductive_bank", "inductive_cable33", "inductive_cable132c"]

"The first `N` loads of `data`, sorted, as harmonic units."
paper_units(data, N = NS) = sort(parse.(Int, collect(keys(data["load"]))))[1:N]

"""
`T` converter types with mean firing delays spread evenly over [0°, 60°), assigned
cyclically over the units: delay 0 is type A (no firing-delay random variable),
every other type has its own random variable shared by its units (30° wide).
WORKING ASSUMPTION for the controlled types, see converter_types.jl.
"""
types_spread(units, T; span = 60.0, w = 30.0) =
    T == 1 ? nothing : types_cyclic(units, [span * (t - 1) / T for t in 1:T]; w = w)

"Replicated data restricted to `units` as harmonic units, local angle reference."
function restricted(data, units)
    hd = HPM.replicate(deepcopy(data), H = H); hd["angle_reference"] = "local"
    for h in harm, k in collect(keys(hd["nw"]["$h"]["load"]))
        parse(Int, k) in units || delete!(hd["nw"]["$h"]["load"], k)
    end
    return hd
end

"PCE parameters for `units`, with optional shared type random variables."
function paper_pce(units, types; bound = "gaussian", eps = PCE_PARAMS["epsilon"])
    p = Dict{String,Any}(deepcopy(PCE_PARAMS)...)
    p["harmonics"] = Dict{Int,Any}(h => rng[h] for h in harm); p["harmonic_units"] = units
    p["bound"] = bound; p["epsilon"] = eps
    types === nothing || (p["type_germs"] = types)
    return p
end

"Hosting current per harmonic of sHHC_SOC for `bound` ∈ cantelli, gaussian, empirical."
function soc_hosting(data, units, types, bound; eps = PCE_PARAMS["epsilon"])
    p = paper_pce(units, types; bound = bound == "empirical" ? "gaussian" : bound, eps = eps)
    sol = bound == "empirical" ?
        HPM.solve_shhc_calibrated(restricted(data, units), sHHC_SOC, mosek; pce_params = p,
                                  hpf_optimizer = ipopt, nsamp = 20000) :
        HPM.solve_shhc(restricted(data, units), sHHC_SOC, mosek; pce_params = p, hpf_optimizer = ipopt)
    return Dict(h => sol["solution"][h] for h in harm)
end

"Harmonic admittance at a (possibly non-integer) harmonic order h, sHHC_SOC stamps."
function ybus(data, h::Real)
    ids, idx = HPM._bus_id_map(data)
    Y = zeros(ComplexF64, length(ids), length(ids))
    for br in values(data["branch"])
        get(br, "br_status", 1) == 0 && continue
        f, t = idx[br["f_bus"]], idx[br["t_bus"]]
        tm = get(br, "tap", 1.0); tm == 0.0 && (tm = 1.0)
        tc = tm * cis(get(br, "shift", 0.0))
        ys = 1 / (br["br_r"] * sqrt(h) + im * br["br_x"] * h)
        yf = get(br, "g_fr", 0.0) + im * get(br, "b_fr", 0.0) * h
        yt = get(br, "g_to", 0.0) + im * get(br, "b_to", 0.0) * h
        Y[f, f] += (ys + yf) / abs2(tc); Y[t, t] += ys + yt
        Y[f, t] -= ys / conj(tc);        Y[t, f] -= ys / tc
    end
    for sh in values(get(data, "shunt", Dict{String,Any}()))
        bs = sh["bs"]
        Y[idx[sh["shunt_bus"]], idx[sh["shunt_bus"]]] += sh["gs"] + im * (bs >= 0 ? bs * h : bs / h)
    end
    for g in values(data["gen"])
        rsc, xsc = get(g, "rsc", 0.0), get(g, "xsc", 0.0)
        (iszero(rsc) && iszero(xsc)) && continue
        Y[idx[g["gen_bus"]], idx[g["gen_bus"]]] += 1 / (rsc * sqrt(h) + im * xsc * h)
    end
    return Y, idx
end

tot(r, k) = sum(r[h][k] for h in harm)
alphas(r) = [r[h][:q95_alpha] for h in harm]
csvline(io, xs...) = (println(io, join(xs, ",")); flush(io))
