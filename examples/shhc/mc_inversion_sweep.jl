################################################################################
# Investigation: MC-vs-sHHC placement for the SOC model on IEEE-30, sweeping the  #
# number of harmonic units N. For each N produce the MC dHHC objective            #
# distribution, the SOC hosting under empirical/Gaussian/Cantelli bounds, and the #
# deterministic dHHC@E[φ] — written as ieee30_N{n}_stoch_* and plotted as          #
# ieee30_N{n}_mc_overlay via plot_mc_overlay.jl.                                   #
#                                                                               #
# Question under test: where does the maximised stochastic hosting fall in the    #
# distribution of per-realisation achievable hosting H(ω)?                         #
################################################################################
using HarmonicPowerModels, PowerModels
using Ipopt, MosekTools, HSL_jll
using JuMP, Printf, Statistics

include(joinpath(@__DIR__, "pce_params.jl"))
include(joinpath(@__DIR__, "shhc_equivalence.jl"))

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

N_LIST = [parse(Int,x) for x in split(get(ENV,"SHHC_NLIST","2,4,8,13,21"), ",")]

summary = NamedTuple[]
for N in N_LIST
    cfg = (; data,H,harm, units=all_units[1:N],
        base_ranges=Dict(h=>(PCE_PARAMS["harmonics"][h].v_min,PCE_PARAMS["harmonics"][h].v_max) for h in harm),
        α̂=PCE_PARAMS["alpha_hat"],β̂=PCE_PARAMS["beta_hat"],ε=PCE_PARAMS["epsilon"],pce_params=PCE_PARAMS,
        ipopt=ipopt_ma27, mosek=mosek_soc, name="ieee30_N$N", outdir=outdir,
        restrict_units=true, qcqp_ok=false, nmc=2000, nsamp_lambda=20000)
    r = stochastic_equivalence(cfg; bounds=["empirical","gaussian","cantelli"])
    o = r.objmc; q(p)=quantile(o,p); vp(h)=100*count(<(h),o)/length(o)
    HE = sum(get_I(solve_dhhc_pm(cfg)[1],h) for h in harm)
    push!(summary, (N=N, q05=q(.05), median=q(.5), q95=q(.95),
                    emp=r.res["empirical"].soc, gau=r.res["gaussian"].soc, can=r.res["cantelli"].soc,
                    dHHC=HE, viol_emp=vp(r.res["empirical"].soc), viol_det=vp(HE)))
end

println("\n", repeat("=",100))
println("SOC sHHC vs MC dHHC distribution — IEEE-30, varying N   (viol = P[H(ω) < hosting] = MC percentile)")
println(repeat("=",100))
@printf("%-4s %-12s %-12s %-12s | %-12s %-9s | %-12s %-9s | %-12s\n",
        "N","MC q05","MC median","MC q95","SOC emp","viol%","dHHC@E[φ]","viol%","SOC cant")
for s in summary
    @printf("%-4d %-12.5e %-12.5e %-12.5e | %-12.5e %-9.1f | %-12.5e %-9.1f | %-12.5e\n",
            s.N,s.q05,s.median,s.q95,s.emp,s.viol_emp,s.dHHC,s.viol_det,s.can)
end
println(repeat("=",100))
println("Reading: viol% = fraction of MC realisations whose max hosting H(ω) is BELOW the committed")
println("hosting (= EN50160 violation probability). Target for the chance constraint is ε=5%.")
