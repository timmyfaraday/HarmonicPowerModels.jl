# diag_vresp.jl — compare LP vs QCQP voltage response C_re[n,k] per PCE index k
# (bus 8, h=11), normalised per unit hosting current. Where do they diverge?
using HarmonicPowerModels, PowerModels
using Ipopt, MosekTools, HSL_jll, JuMP, Printf
const PMs=PowerModels; const HPM=HarmonicPowerModels
include(joinpath(HPM.BASE_DIR,"examples","shhc","pce_params.jl"))
ipopt=optimizer_with_attributes(Ipopt.Optimizer,"linear_solver"=>"ma27","hsllib"=>HSL_jll.libhsl_path,"print_level"=>0)
mosek=optimizer_with_attributes(Mosek.Optimizer,"QUIET"=>true)
data=PMs.parse_file(joinpath(HPM.BASE_DIR,"examples","shhc","data","industrial_network_shhc.m"))
H=[1,5,7,11,13]

hl=HPM.replicate(deepcopy(data),H=H)
sl=HPM.solve_shhc(hl, sHHC_LP, mosek; pce_params=deepcopy(PCE_PARAMS), hpf_optimizer=ipopt)
pce=hl["pce"]; aRe=pce["alpha_re"]; aIm=pce["alpha_im"]; ind=pce["basis"].ind
fund=hl["nw"]["1"]; bus_ids=sort([b["bus_i"] for (_,b) in fund["bus"]]); bidx=Dict(id=>i for (i,id) in enumerate(bus_ids))

hq=HPM.replicate(deepcopy(data),H=H)
sq=HPM.solve_shhc(hq, sHHC_QCQP, ipopt; pce_params=deepcopy(PCE_PARAMS))
Ih = sq["solution"]["nw"]["11"]["load"]["1"]["cmd"]
vr = sq["solution"]["nw"]["11"]["bus"]["8"]["vr_pce"]   # Dict k=>vr_pce
vi = sq["solution"]["nw"]["11"]["bus"]["8"]["vi_pce"]

n8=bidx[8]; P=pce["basis"].dim
println("h=11 bus 8: voltage response Cre[k] (LP, per-I) vs vr_pce[k]/Ih (QCQP)")
@printf("%-4s %-10s %-14s %-14s %-14s %-14s\n","k","multi-idx","Cre_LP","Cre_QCQP","Cim_LP","Cim_QCQP")
for k in 0:P-1
    cl=aRe[11][n8,k+1]; il=aIm[11][n8,k+1]
    cq=get(vr,string(k),get(vr,k,NaN))/Ih; iq=get(vi,string(k),get(vi,k,NaN))/Ih
    (abs(cl)>1e-12 || abs(cq)>1e-12) || continue
    mi = join(ind[k+1,:],",")
    @printf("%-4d (%-7s) %-14.6e %-14.6e %-14.6e %-14.6e\n",k,mi,cl,cq,il,iq)
end
