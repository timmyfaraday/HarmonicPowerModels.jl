################################################################################
# Does the 6-ideal-ground model suppress recombination vs a single ground with    #
# generator subtransient impedance (the IEEE-30 "991" config)?                    #
#                                                                               #
# Builds the Kron-reduced harmonic transfer impedance Z_h for both grounding      #
# models and, at the most-exposed bus, Monte-Carlo samples the per-unit IHD       #
# response a_h(ω)=|Σ_d Z[n,d]e^{jφ_d}|². The recombination potential is the gap   #
# between the coherent worst case and the q95 of a_h: R = a_worst / q95(a) (the    #
# factor by which a 95%-reliable hosting beats the robust worst case). Bigger R    #
# ⇒ more recombination headroom. Independent germ per unit; section-II ranges.     #
################################################################################
using HarmonicPowerModels, PowerModels
using LinearAlgebra, Distributions, Random, Statistics, Printf
const HPM = HarmonicPowerModels
include(joinpath(@__DIR__, "pce_params.jl"))

data = PowerModels.parse_file(joinpath(@__DIR__, "data", "30bus_network_shhc.m"))
H = [1,5,7,11,13]; harm = [5,7,11,13]
α̂, β̂ = PCE_PARAMS["alpha_hat"], PCE_PARAMS["beta_hat"]; ε = PCE_PARAMS["epsilon"]
base = Dict(h => (PCE_PARAMS["harmonics"][h].v_min, PCE_PARAMS["harmonics"][h].v_max) for h in harm)
hd = HPM.replicate(deepcopy(data), H = H)

genbuses = [1,2,5,8,11,13]
# generator subtransient impedance (from the data-file header), fundamental values
gen_z = Dict(1=>(0.0025,0.05), 2=>(0.0025,0.05), 5=>(0.025,0.5),
             8=>(0.025,0.5), 11=>(0.025,0.5), 13=>(0.025,0.5))

busids = sort([b["bus_i"] for (_,b) in hd["nw"]["1"]["bus"]])
bidx   = Dict(id=>i for (i,id) in enumerate(busids)); Nn = length(busids)
units  = sort(parse.(Int, collect(keys(hd["nw"]["1"]["load"]))))           # all 21
ubus   = Dict(d => hd["nw"]["1"]["load"]["$d"]["load_bus"] for d in units)

# branch stamp (tap-aware) into Y for harmonic network ntw
function stamp_branches!(Y, ntw)
    for (_,br) in ntw["branch"]
        f=bidx[br["f_bus"]]; t=bidx[br["t_bus"]]; r=br["br_r"]; x=br["br_x"]
        (abs(r)+abs(x)>0)||continue
        ys=1/(r+im*x); tm=get(br,"tap",1.0); tm==0 && (tm=1.0); sh=get(br,"shift",0.0)
        tc=tm*cis(sh); bf=get(br,"b_fr",0.0); bt=get(br,"b_to",0.0)
        Y[f,f]+=(ys+im*bf)/abs2(tc); Y[t,t]+=ys+im*bt
        Y[f,t]-=ys/conj(tc);          Y[t,f]-=ys/tc
    end
end

# reduced Z for harmonic h under a grounding model; returns (Z, pos) over non-ground buses
function Zmodel(h; single_ground::Bool)
    ntw = hd["nw"]["$h"]
    if !single_ground
        Y = zeros(ComplexF64, Nn, Nn); stamp_branches!(Y, ntw)
        gnd = Set(bidx[g] for g in genbuses)                 # 6 ideal grounds
        nonref = [i for i in 1:Nn if i ∉ gnd]
    else
        # add fictitious ground node 991; gens connect via X'' (scaled: r√h, x·h)
        M = Nn + 1; g991 = M
        Y = zeros(ComplexF64, M, M); stamp_branches!(Y, ntw)
        for g in genbuses
            r,x = gen_z[g]; rh = r*sqrt(h); xh = x*h; ys = 1/(rh+im*xh)
            f = bidx[g]
            Y[f,f]+=ys; Y[g991,g991]+=ys; Y[f,g991]-=ys; Y[g991,f]-=ys
        end
        gnd = Set([g991])                                    # single ground
        nonref = [i for i in 1:M if i ∉ gnd]
    end
    Yr = Y[nonref, nonref]
    Z  = Yr \ Matrix{ComplexF64}(I, length(nonref), length(nonref))
    pos = Dict(busids[i]=>k for (k,i) in enumerate(nonref) if i<=Nn)
    return Z, pos
end

Random.seed!(1); Bd = Distributions.Beta(α̂,β̂); nsamp = 40_000
λstar(x)=(σ=std(x); σ>0 ? (quantile(x,1-ε)-mean(x))/σ : 0.0)

function analyse(single_ground)
    # sample germ per unit once, reuse across harmonics (comonotone)
    Ω = [Dict(d=>rand(Bd) for d in units) for _ in 1:nsamp]
    out = Dict{Int,NamedTuple}()
    for h in harm
        Z,pos = Zmodel(h; single_ground=single_ground)
        # most-exposed non-ground bus = max Σ_d |Z[n,d]|
        exposure(n) = sum(haskey(pos,ubus[d]) ? abs(Z[pos[n],pos[ubus[d]]]) : 0.0 for d in units)
        nstar = argmax(Dict(n=>exposure(n) for n in keys(pos)))
        zc = Dict(d=>(haskey(pos,ubus[d]) ? Z[pos[nstar],pos[ubus[d]]] : 0im) for d in units)
        a_wc = (sum(abs(zc[d]) for d in units))^2            # coherent worst |V|²/I²
        a = Vector{Float64}(undef,nsamp)
        lo=base[h][1]; wid=base[h][2]-base[h][1]
        for m in 1:nsamp
            Vr=0.0;Vi=0.0
            for d in units
                φ=lo+wid*Ω[m][d]; zr,zi=real(zc[d]),imag(zc[d]); c,s=cos(φ),sin(φ)
                Vr+=zr*c-zi*s; Vi+=zi*c+zr*s
            end
            a[m]=Vr^2+Vi^2
        end
        out[h]=(nstar=nstar, a_wc=a_wc, mean=mean(a), q05=quantile(a,.05), q95=quantile(a,.95),
                R=a_wc/quantile(a,.95), λ=λstar(a))
    end
    out
end

println("="^96)
println("Recombination potential by grounding model — IEEE-30, all $(length(units)) units, n=$nsamp")
println("R = a_worst/q95(a): hosting headroom of chance vs robust worst case (∝ recombination).")
println("="^96)
for (tag, sg) in (("6 ideal grounds (current)", false), ("1 ground + gen X'' (991)", true))
    r = analyse(sg)
    println("\n$tag:")
    @printf("  %-4s %-6s %-12s %-12s %-12s %-8s %-8s\n","h","bus*","a_worst","E[a]","q95(a)","R","λ*")
    for h in harm
        x=r[h]; @printf("  %-4d %-6d %-12.4e %-12.4e %-12.4e %-8.3f %-8.3f\n",
                        h, x.nstar, x.a_wc, x.mean, x.q95, x.R, x.λ)
    end
end
println("\n", "="^96)
println("Higher R (and lower λ*) under the single-ground model ⇒ generators-as-impedance unlocks")
println("more recombination than 6 ideal grounds.")
