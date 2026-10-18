# HarmonicPowerModels.jl                                                       #
# Phase-2 validation: compare the rewired (PolyChaos / PCEBasis) LP-reduction  #
# pipeline against the Phase-0 golden reference.                               #
#                                                                               #
# Only basis-INDEPENDENT quantities are compared (I_star, binding buses,       #
# C_scalar vectors): both implementations project onto the same polynomial     #
# subspace w.r.t. the same probability measure, so these must agree to         #
# rounding. Raw coefficient vectors are basis-dependent and are NOT compared.  #
#                                                                               #
# Run: julia --project test/golden/compare_phase2_golden.jl                    #
# Author: Tom Van Acker                                                         #

using HarmonicPowerModels, PowerModels
using Printf, Test

include(joinpath(@__DIR__, "..", "..", "examples", "shhc", "pce_params.jl"))

const PMs  = PowerModels
const HPM  = HarmonicPowerModels
const HARM = [5, 7, 11, 13]

"Parse `name = v1 v2 …` lines of a golden file into a Dict."
function parse_golden(path)
    d = Dict{String, Vector{Float64}}()
    for line in eachline(path)
        (startswith(line, "#") || !occursin(" = ", line)) && continue
        name, vals = split(line, " = "; limit = 2)
        d[name] = parse.(Float64, split(vals))
    end
    return d
end

function compare(tag::String, casefile::String)
    g    = parse_golden(joinpath(@__DIR__, "pce_golden_$(tag).txt"))
    data = PMs.parse_file(joinpath(@__DIR__, "..", "..", "examples", "shhc",
                                   "data", casefile))

    beta_params = (alpha_shape = PCE_PARAMS["alpha_hat"],
                   beta_shape  = PCE_PARAMS["beta_hat"],
                   delta       = PCE_PARAMS["delta"])
    ε = PCE_PARAMS["epsilon"]

    unit_ids = sort(parse.(Int, collect(keys(data["load"]))))

    @testset "$tag vs Phase-0 golden" begin
        for h in HARM
            hb = PCE_PARAMS["harmonics"][h]
            angle_params_h = Dict(uid => (rad2deg(hb.v_min), rad2deg(hb.v_max))
                                  for uid in unit_ids)
            r = HPM.compare_shared_vs_independent(data, h, ε, beta_params,
                                                  angle_params_h)

            @test r.I_star_independent ≈ g["h$(h).I_star_independent"][1] rtol = 1e-9
            @test r.I_star_shared      ≈ g["h$(h).I_star_shared"][1]      rtol = 1e-9
            @test r.binding_bus_independent == Int(g["h$(h).binding_bus_independent"][1])
            @test r.binding_bus_shared      == Int(g["h$(h).binding_bus_shared"][1])
            @test r.C_scalar_independent ≈ g["h$(h).C_scalar_independent"] rtol = 1e-9
            @test r.C_scalar_shared      ≈ g["h$(h).C_scalar_shared"]      rtol = 1e-9

            @printf("  %s h=%-2d  I_ind: %.10e vs %.10e   I_sh: %.10e vs %.10e\n",
                    tag, h, r.I_star_independent, g["h$(h).I_star_independent"][1],
                    r.I_star_shared, g["h$(h).I_star_shared"][1])
        end
    end
end

compare("industrial", "industrial_network_shhc.m")
compare("ieee30",     "30bus_network_shhc.m")
println("\n✓ Phase-2 pipeline matches Phase-0 goldens on all basis-independent quantities.")
