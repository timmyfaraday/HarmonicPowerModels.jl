################################################################################
# Pinpoint the dHHC↔reduced-Z gap: compare dHHC native voltages to V=Z·I where    #
# Z is built (a) WITH branch pi-shunt (as _harmonic_admittance does) and          #
# (b) WITHOUT shunt — using dHHC's own currents. Whichever matches isolates the    #
# discrepancy.                                                                     #
################################################################################
using HarmonicPowerModels, PowerModels, Ipopt, MosekTools, HSL_jll, JuMP, Printf, LinearAlgebra
const HPM = HarmonicPowerModels
include(joinpath(@__DIR__,"pce_params.jl"))
ipopt=optimizer_with_attributes(Ipopt.Optimizer,"linear_solver"=>"ma27","hsllib"=>HSL_jll.libhsl_path,"print_level"=>0)
mosek=optimizer_with_attributes(Mosek.Optimizer,"QUIET"=>true)
data=PowerModels.parse_file(joinpath(@__DIR__,"data","30bus_network_shhc.m"))
H=[1,5,7,11,13]; harm=[5,7,11,13]
Eω=PCE_PARAMS["alpha_hat"]/(PCE_PARAMS["alpha_hat"]+PCE_PARAMS["beta_hat"])
Eφd(h)=(hb=PCE_PARAMS["harmonics"][h]; rad2deg(hb.v_min+(hb.v_max-hb.v_min)*Eω))
units=sort(parse.(Int,collect(keys(data["load"]))))[1:4]; uset=Set(units)
hd=HPM.replicate(deepcopy(data),H=H)
for h in harm,k in collect(keys(hd["nw"]["$h"]["load"])); parse(Int,k) in uset ? (hd["nw"]["$h"]["load"][k]["ref_angle"]=Eφd(h)) : delete!(hd["nw"]["$h"]["load"],k); end
sol=HPM.solve_hhc(hd,dHHC_SOC,mosek,ipopt)

fund=hd["nw"]["1"]; busids=sort([b["bus_i"] for (_,b) in fund["bus"]])
bidx=Dict(id=>i for (i,id) in enumerate(busids)); Nn=length(busids)
refset=Set(bidx[b["bus_i"]] for (_,b) in fund["bus"] if b["bus_type"]==3)
nonref=[i for i in 1:Nn if i∉refset]; nrid=[busids[i] for i in nonref]
bpos=Dict(id=>k for (k,id) in enumerate(nrid))
ubus=Dict(d=>fund["load"]["$d"]["load_bus"] for d in units)

function Zof(nw; shunt::Bool, tap::Bool=false)
    Y=zeros(ComplexF64,Nn,Nn)
    for (_,br) in nw["branch"]
        f=bidx[br["f_bus"]]; t=bidx[br["t_bus"]]; r=br["br_r"]; x=br["br_x"]
        (abs(r)+abs(x)>0) || continue
        ys=1/(r+im*x)
        tm = tap ? get(br,"tap",1.0) : 1.0
        sh = tap ? get(br,"shift",0.0) : 0.0
        tc = tm*cis(sh)                                   # complex tap t = tm·e^{jθ}
        bf=shunt ? im*get(br,"b_fr",0.0) : 0.0im
        bt=shunt ? im*get(br,"b_to",0.0) : 0.0im
        Y[f,f]+=(ys+bf)/abs2(tc); Y[t,t]+=ys+bt
        Y[f,t]-=ys/conj(tc);      Y[t,f]-=ys/tc
    end
    Yr=Y[nonref,nonref]; return Yr\Matrix{ComplexF64}(I,length(nonref),length(nonref))
end

function meandev(Z)
    tot=Float64[]
    for h in harm
        I=Dict(d=>complex(sol["solution"]["nw"]["$h"]["load"]["$d"]["crd"],
                          sol["solution"]["nw"]["$h"]["load"]["$d"]["cid"]) for d in units)
        Zh=Z(hd["nw"]["$h"])
        for (bi,b) in sol["solution"]["nw"]["$h"]["bus"]
            n=parse(Int,bi); haskey(bpos,n)||continue
            Vd=complex(b["vr"],b["vi"]); abs(Vd)<1e-9 && continue
            Vz=sum((haskey(bpos,ubus[d]) ? Zh[bpos[n],bpos[ubus[d]]] : 0.0im)*I[d] for d in units)
            push!(tot, abs(abs(Vz)-abs(Vd))/abs(Vd)*100)
        end
    end
    sum(tot)/length(tot)
end

@printf("mean |Δ|V|| shunt, no tap   : %.4f %%\n", meandev(nw->Zof(nw;shunt=true,tap=false)))
@printf("mean |Δ|V|| shunt + TAP     : %.4f %%\n", meandev(nw->Zof(nw;shunt=true,tap=true)))
