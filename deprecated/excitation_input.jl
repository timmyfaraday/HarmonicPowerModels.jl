## BH-CURVE EXAMPLE
# pkgs
using Dierckx

# Thyssenkrupp - PowerCore H 100-23 50Hz
B = [0.000, 0.144, 0.200, 0.260, 0.328, 0.400, 0.504, 0.600, 0.695, 1.528, 1.716, 1.776, 1.816, 1.828, 1.832, 1.845, 1.856, 1.860]
H = [0.000, 3.000, 4.000, 5.000, 6.000, 7.000, 8.000, 9.000, 10.00, 20.00, 30.00, 40.00, 50.00, 60.00, 70.00, 80.00, 90.00, 100.0]
BH_powercore_h100_23 = Dierckx.Spline1D(B, H; k=3, bc="nearest")

# Input
xfmr_exc = Dict(1 => Dict(  "Hᴱ"    => [1,3],
                            "Hᴵ"    => [1,3],
                            "Fᴱ"    => :rectangular,
                            "Fᴵ"    => :rectangular,
                            "l"     => 0.2,
                            "A"     => 1.0,
                            "N"     => 1000,
                            "BH"    => BH_powercore_h100_23)
                )
