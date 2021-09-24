# 2-bus
freq = 50 #Hz
harmonics = 1:2:3 #h in {1,3}

function thd(V)
    return hypot(abs.(V[2:end]...))/abs(V[1])
end

function complexpower(V,I)
    return V*conj(I)
end

# reference bus voltage (i=0)
V0 = zeros(length(harmonics))
V0[1] = 1 + 0*im # 1 V at 50 Hz, no voltage except for fundamental frequency

# branch as a series impedance
r_branch_h1 = 1 #ohm
x_branch_h1 = 1 #ohm
r_branch_h3 = r_branch_h1 * sqrt(3) #ohm
x_branch_h3 = 3 * x_branch_h1 #ohm 
z_branch = [r_branch_h1 + im*x_branch_h1; r_branch_h3 + im*x_branch_h3]

# constant current load on bus i=1
i_load_h1 = 0.1 - 0.1*im;
i_load_h3 = 0.01 + 0*im;
i_load = [ i_load_h1; i_load_h3]

## load current = current provided by generator on bus 0
i_01 = i_load
s_01 = complexpower.(V0,i_01) 

# complex power flowing in the branch in the direction of bus 0 to bus 1
@show s_01

# solve Ohm's law to obtain voltage at bus 1
V1 = V0 .- z_branch.* i_01
s_10 = complexpower.(V1,-i_01) 

# complex power flowing in the branch in the direction of bus 1 to bus 0
@show s_10

#complex losses in branch
s_01_loss = s_01 + s_10
#complex losses at fundamental
@show s_01_loss[1]
#complex losses due to harmonics
s_01_loss_harmonics = sum(real(s_01_loss[2:end]))
@show s_01_loss_harmonics


# power consumed by constant current load
s_load = complexpower.(V1,i_load) 
@show s_load

# voltage magnitude of the harmonics
@show abs.(V1)

# total harmonic distortion
@show thd(V1)

