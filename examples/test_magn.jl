using Plots
using Dierckx
using SignalDecomposition
using HarmonicPowerModels
using PowerModels
using Test

Sbase = 100e6
Vbase = 12470
H     = [1,3]
Hᴵ    = collect(1:2:100)
Ere   = [1.0,0.0]
Eim   = [0.0,0.0]
l     = 10.0
A     = 0.5
N     = 75

w     = 2 .* pi .* 50.0 .* H
t     = 0.0:0.0001:0.1

Abase = Sbase / Vbase

B⁺ = [0.144, 0.200, 0.260, 0.328, 0.400, 0.504, 0.600, 0.695, 1.528, 1.716, 1.776, 1.816, 1.828, 1.832, 1.845, 1.856, 1.860]
H⁺ = [3.000, 4.000, 5.000, 6.000, 7.000, 8.000, 9.000, 10.00, 20.00, 30.00, 40.00, 50.00, 60.00, 70.00, 80.00, 90.00, 100.0]
BH_powercore_h100_23 = Dierckx.Spline1D(vcat(reverse(-B⁺),0.0,B⁺), vcat(reverse(-H⁺),0.0,H⁺); k=3, bc="nearest")

# e(t) = |E| ⋅ sin(w t + ∠E)
e   = sqrt(2) .* Vbase .* sum(hypot(Ere[ni],Eim[ni]) .* sin.(w[ni] .* t .+ atan(Eim[ni],Ere[ni])) for ni in 1:length(H))

plot(t,e)

# b(t) = 1 / AN ∫ -e(t) dt = |E| / (w A N) * - (-cos(w t + ∠E)) = |E| / (w A N) * cos(w t + ∠E)
b   = Vbase ./ (A .* N).* sum(hypot(Ere[ni],Eim[ni]) ./ w[ni] .* cos.(w[ni] .* t .+ atan(Eim[ni],Ere[ni])) for ni in 1:length(H))

plot(t,b)

# h(t) = BH(b(t))
h   = BH_powercore_h100_23.(b)

plot(t,h)

# i = l / N * h
i   = l ./ N .* h 

plot(t,i)

I   = i ./ Abase 

plot(t,I)

fq = Sinusoidal(50.0 .* Hᴵ)
decompose(collect(t), I, fq)

IM = sum(fq.A[ni] .* cos.(2.0 .* pi .* 50.0 .* t .+ fq.φ[ni]) for ni in 2:length(fq.A))

## we take sin as basis, therefore the angle needs to be shifted with 90 degrees to be correct!

plot!(t,IM)

## Alternatief
# pkg const
const PMs = PowerModels
const HPM = HarmonicPowerModels

# path to the data
path = joinpath(HPM.BASE_DIR,"test/data/matpower/xfmr/case_xfmr_YNyn0.m")

# load data
data = PMs.parse_file(path)
# xfmr magnetizing data
magn = Dict("Hᴱ"    => [1,3], 
            "Hᴵ"    => collect(1:2:25),
            "Fᴱ"    => :rectangular,
            "Fᴵ"    => :rectangular,
            "xfmr"  => Dict(1 => Dict(  "l"  => 10.0,
                                        "A"  => 0.5,
                                        "N"  => 75,
                                        "BH" => BH_powercore_h100_23)
                            )
            )

# harmonic data
hdata = HPM.replicate(data, xfmr_magn=magn)

# get first harmonic amplitude through first approach
A1_first = fq.A[2] / sqrt(2)
Φ1_first = fq.φ[2]

# get first harmonic amplitude through interpolate approach
xfmr = hdata["nw"]["1"]["xfmr"]["1"]
A1_second = hypot(xfmr["Im_A"].(1.0,0.0,0.0,0.0), xfmr["Im_B"].(Ere...,Eim...))
Φ1_second = atan(xfmr["Im_B"].(Ere...,Eim...), xfmr["Im_A"].(Ere...,Eim...))

@test isapprox(A1_first, A1_second, atol=1e-5)

