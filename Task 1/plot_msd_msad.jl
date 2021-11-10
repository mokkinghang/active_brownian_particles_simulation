include("Task1.jl")

N1, N2 = 40,25
L = 85
dt = 10^(-3)
writedata_steps = 10
T = 100
sigma_to_R = 1
enable_pbc = false

N = N1*N2
time_steps = trunc(Int,T/dt)
M = trunc(Int,T/(dt*writedata_steps)) + 1

Pe_arr = [0,5,10,15,20]

save_traj_to_jld2(Pe_arr, enable_pbc=enable_pbc)

plot_msd(["Pe=0.jld2", "Pe=5.jld2", "Pe=10.jld2", "Pe=15.jld2", "Pe=20.jld2"])
plot_msad(["Pe=0.jld2", "Pe=5.jld2", "Pe=10.jld2", "Pe=15.jld2", "Pe=20.jld2"])
