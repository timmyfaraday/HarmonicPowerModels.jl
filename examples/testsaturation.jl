
using SignalDecomposition
using Interpolations
using Plots

const _SDC = SignalDecomposition
const _INT = Interpolations

const freq = 50.0
const nw_id_default = 1

excitation_flux_polar(V, θ, w, t) = 
    sum(V[n] ./ w[n] .* sin.(w[n] .* t .+ θ[n]) for n in 1:length(V))
excitation_flux_rectangular(Vre, Vim, w, t) =
    excitation_flux_polar(hypot.(Vre,Vim), atan.(Vim,Vre), w, t)
excitation_current_sigmoid(inom, ψmax, ψ) = 
    -inom .* log.(2.0 ./ (ψ ./ ψmax .+ 1.0) .- 1.0)

""
function sample_voltage_polar(voltage_harmonics, dv, vmin, vmax, dθ, θmin, θmax)
    S = reduce(vcat,[[vmin[ni]:dv[ni]:vmax[ni],θmin[ni]:dθ[ni]:θmax[ni]] 
                      for (ni,nh) in enumerate(voltage_harmonics)])
    R = [1:length(s) for s in S]
    return S, R
end
function sample_voltage_rectangular(voltage_harmonics, dv, vmin, vmax)
    S = reduce(vcat,[[vmin[ni]:dv[ni]:vmax[ni],vmin[ni]:dv[ni]:vmax[ni]] 
                      for (ni,nh) in enumerate(voltage_harmonics)])
    R = [1:length(s) for s in S]
    return S, R
end

xfmr_exc = xfmr

## HARMONICS
voltage_harmonics = xfmr_exc["voltage_harmonics"]
current_harmonics = xfmr_exc["current_harmonics"]

## TIME
N = xfmr_exc["N"]

## EXCITATION CURRENT 
current_type = xfmr_exc["current_type"]

## EXCITATION FLUX 
excitation_type = xfmr_exc["excitation_type"] 
inom = xfmr_exc["inom"]
ψmax = xfmr_exc["ψmax"]

## EXCITATION VOLTAGE 
voltage_type = xfmr_exc["voltage_type"]
dv   = xfmr_exc["dv"]
vmin = xfmr_exc["vmin"]
vmax = xfmr_exc["vmax"]
dθ   = xfmr_exc["dθ"]
θmin = xfmr_exc["θmin"]
θmax = xfmr_exc["θmax"]

## TIME
w = (2.0 * pi * freq) .* voltage_harmonics
t = 0.0:(1 / (freq * N * maximum(current_harmonics))):(5.0 / freq)

## DECOMPOSITION 
fq = _SDC.Sinusoidal(freq .* current_harmonics)

## SAMPLE VOLTAGE
if voltage_type == :polar
    S, R = sample_voltage_polar(voltage_harmonics, dv, vmin, vmax, dθ, θmin, θmax)
elseif voltage_type == :rectangular
    S, R = sample_voltage_rectangular(voltage_harmonics, dv, vmin, vmax)
end
Ia = Dict(nh => zeros(R...) for nh in current_harmonics)
Ib = Dict(nh => zeros(R...) for nh in current_harmonics)

nr = [20,20,12,12]
sample = [S[ni][ns] for (ni,ns) in enumerate(nr)]

# if voltage_type == :polar
#     V, θ = sample[1:2:end], sample[2:2:end]
#     ψexc = excitation_flux_polar(V, θ, w, t)
# else voltage_type == :rectangular 
Vre, Vim = sample[1:2:end], sample[2:2:end]
ψexc = excitation_flux_rectangular(Vre, Vim, w, t) 
# end

# if excitation_type == :sigmoid
I_exc = excitation_current_sigmoid(inom, ψmax, ψexc)
# end


plot(t,ψexc)
plot!(t,I_exc)

_SDC.decompose(t, I_exc, fq)

I, φ = fq.A[2:end], -fq.φ[2:end]