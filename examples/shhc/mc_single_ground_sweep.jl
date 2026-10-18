################################################################################
# Rerun the MC/SOC inversion sweep on the SINGLE-GROUND model (bus 1 fundamental  #
# slack kept; generators reach harmonic ground only via their X'' impedance) and  #
# compare against the current 6-ideal-ground model. Produces ieee30_1g_N{n}_*     #
# outputs + mc_overlay figures (deterministic dHHC@E[φ] marker included).          #
################################################################################
using HarmonicPowerModels, PowerModels
using Ipopt, MosekTools, HSL_jll
using JuMP, Printf, Statistics

include(joinpath(@__DIR__, "pce_params.jl"))
include(joinpath(@__DIR__, "shhc_equivalence.jl"))
include(joinpath(@__DIR__, "single_ground.jl"))

const TLIM = 360.0
ipopt_ma27 = optimizer_with_attributes(Ipopt.Optimizer,
                "linear_solver"=>"ma27","hsllib"=>HSL_jll.libhsl_path,"print_level"=>0,
                "max_wall_time"=>TLIM,"tol"=>1e-10,"constr_viol_tol"=>1e-10,
                "acceptable_tol"=>1e-10,"acceptable_constr_viol_tol"=>1e-10)
mosek_soc  = optimizer_with_attributes(Mosek.Optimizer,"QUIET"=>true,
                "MSK_DPAR_OPTIMIZER_MAX_TIME"=>TLIM)

data = PowerModels.parse_file(joinpath(@__DIR__,"data","30bus_network_shhc.m"))
H=[1,5,7,11,13]; harm=[5,7,11,13]
outdir = joinpath(@__DIR__,"output"); mkpath(outdir)
all_units = sort(parse.(Int, collect(keys(data["load"]))))
N_LIST = [parse(Int,x) for x in split(get(ENV,"SHHC_NLIST","4,8,13,21"), ",")]

summary = NamedTuple[]
for N in N_LIST
    cfg = (; data,H,harm, units=all_units[1:N],
        base_ranges=Dict(h=>(PCE_PARAMS["harmonics"][h].v_min,PCE_PARAMS["harmonics"][h].v_max) for h in harm),
        α̂=PCE_PARAMS["alpha_hat"],β̂=PCE_PARAMS["beta_hat"],ε=PCE_PARAMS["epsilon"],pce_params=PCE_PARAMS,
        ipopt=ipopt_ma27, mosek=mosek_soc, name="ieee30_1g_N$N", outdir=outdir,
        restrict_units=true, qcqp_ok=false, nmc=2000, nsamp_lambda=20000,
        postrep = to_single_ground!)
    r = stochastic_equivalence(cfg; bounds=["empirical","gaussian","cantelli"])
    o = r.objmc; q(p)=quantile(o,p); vp(h)=100*count(<(h),o)/length(o)
    HE = sum(get_I(solve_dhhc_pm(cfg)[1],h) for h in harm)
    push!(summary, (N=N, q05=q(.05), median=q(.5), q95=q(.95), spread=q(.95)/q(.05),
                    emp=r.res["empirical"].soc, dHHC=HE, viol_emp=vp(r.res["empirical"].soc)))
end

println("\n", repeat("=",96))
println("SINGLE-GROUND (gen X'') model — SOC sHHC vs MC dHHC, varying N   [compare to 6-ground sweep]")
println(repeat("=",96))
@printf("%-4s %-12s %-12s %-12s %-9s | %-12s %-9s | %-12s\n",
        "N","MC q05","MC median","MC q95","q95/q05","SOC emp","viol%","dHHC@E[φ]")
for s in summary
    @printf("%-4d %-12.5e %-12.5e %-12.5e %-9.3f | %-12.5e %-9.1f | %-12.5e\n",
            s.N,s.q05,s.median,s.q95,s.spread,s.emp,s.viol_emp,s.dHHC)
end
println(repeat("=",96))
println("Wider q95/q05 spread than the 6-ground model ⇒ generators-as-impedance unlocks recombination.")
