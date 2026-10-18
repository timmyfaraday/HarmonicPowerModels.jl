################################################################################
# Diagnostic: locate the ~1% dHHC↔reduced-Z model gap on IEEE-30.                 #
# Take dHHC's own optimal current injections at E[φ] and reconstruct the bus      #
# voltages through the reduced harmonic impedance (V = Z·I). Compare against       #
# dHHC's native voltages per bus/harmonic — divergence ⇒ the reduced-Z network     #
# differs from the full IVR-dHHC network there.                                    #
################################################################################
using HarmonicPowerModels, PowerModels
using Ipopt, MosekTools, HSL_jll, JuMP, Printf
const HPM = HarmonicPowerModels
include(joinpath(@__DIR__, "pce_params.jl"))

ipopt = optimizer_with_attributes(Ipopt.Optimizer,"linear_solver"=>"ma27",
            "hsllib"=>HSL_jll.libhsl_path,"print_level"=>0)
mosek = optimizer_with_attributes(Mosek.Optimizer,"QUIET"=>true)
data = PowerModels.parse_file(joinpath(@__DIR__,"data","30bus_network_shhc.m"))
H=[1,5,7,11,13]; harm=[5,7,11,13]
Eω=PCE_PARAMS["alpha_hat"]/(PCE_PARAMS["alpha_hat"]+PCE_PARAMS["beta_hat"])
Eφd(h)=(hb=PCE_PARAMS["harmonics"][h]; rad2deg(hb.v_min+(hb.v_max-hb.v_min)*Eω))
all_units=sort(parse.(Int,collect(keys(data["load"])))); units=all_units[1:4]; uset=Set(units)

hd=HPM.replicate(deepcopy(data),H=H)
for h in harm, k in collect(keys(hd["nw"]["$h"]["load"]))
    parse(Int,k) in uset ? (hd["nw"]["$h"]["load"][k]["ref_angle"]=Eφd(h)) : delete!(hd["nw"]["$h"]["load"],k)
end
sol=HPM.solve_hhc(hd,dHHC_SOC,mosek,ipopt)

Z,bp=HPM.reduced_harmonic_impedance(hd)
fund=hd["nw"]["1"]
ubus=Dict(d=>fund["load"]["$d"]["load_bus"] for d in units)
refset=Set(b["bus_i"] for (_,b) in fund["bus"] if b["bus_type"]==3)

println("ref (gen) buses: ", sort(collect(refset)))
println("unit buses: ", [ubus[d] for d in units], "  (ref? ", [ubus[d] in refset for d in units], ")")

function compare()
    @printf("\n%-4s %-4s %-13s %-13s %-9s | injection at non-ref unit buses only\n","bus","h","|V|dHHC","|V|=Z·I","Δ%%")
    maxd=0.0
    for h in harm
        I=Dict(d=>complex(sol["solution"]["nw"]["$h"]["load"]["$d"]["crd"],
                          sol["solution"]["nw"]["$h"]["load"]["$d"]["cid"]) for d in units)
        for (bi,b) in sol["solution"]["nw"]["$h"]["bus"]
            n=parse(Int,bi); haskey(bp,n)||continue
            Vd=complex(b["vr"],b["vi"])
            Vz=sum((haskey(bp,ubus[d]) ? Z[h][bp[n],bp[ubus[d]]] : 0.0im)*I[d] for d in units)
            ad=abs(Vd); az=abs(Vz)
            (ad<1e-9 && az<1e-9) && continue
            dd=abs(Vz-Vd)/max(ad,1e-12)*100; maxd=max(maxd,dd)
            dd>0.5 && @printf("%-4d %-4d %-13.5e %-13.5e %-9.3f\n",n,h,ad,az,dd)
        end
    end
    @printf("\nMax |V| deviation (Z·I vs dHHC native, dHHC's own currents): %.3f %%\n", maxd)
end
compare()

# magnitude-only deviation per harmonic (ignores the 180° sign convention)
function permag()
    println("\nMagnitude-only |V| deviation per harmonic (mean over non-ref buses):")
    for h in harm
        I=Dict(d=>complex(sol["solution"]["nw"]["$h"]["load"]["$d"]["crd"],
                          sol["solution"]["nw"]["$h"]["load"]["$d"]["cid"]) for d in units)
        ds=Float64[]
        for (bi,b) in sol["solution"]["nw"]["$h"]["bus"]
            n=parse(Int,bi); haskey(bp,n)||continue
            Vd=complex(b["vr"],b["vi"])
            Vz=sum((haskey(bp,ubus[d]) ? Z[h][bp[n],bp[ubus[d]]] : 0.0im)*I[d] for d in units)
            abs(Vd)<1e-9 && continue
            push!(ds, abs(abs(Vz)-abs(Vd))/abs(Vd)*100)
        end
        @printf("  h=%-3d mean |Δ|V|| = %.3f %%  (n=%d buses)\n", h, sum(ds)/length(ds), length(ds))
    end
end
permag()

# branch series r,x and shunt b at each harmonic (frequency scaling check)
println("\nSample branch series/shunt at h=5 vs h=13 (scaled in replicate):")
for h in (5,13)
    nb=0; tb=0.0
    for (_,br) in hd["nw"]["$h"]["branch"]
        b=get(br,"b_fr",0.0); tb+=abs(b); b!=0 && (nb+=1)
    end
    @printf("  h=%-3d: branches with b_fr≠0 = %d,  Σ|b_fr| = %.4e\n", h, nb, tb)
end
println("(small ⇒ reduced-Z matches dHHC network; large ⇒ genuine model difference)")
