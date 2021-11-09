include("Task1.jl")

N1, N2 = 1,5
L = 85
dt = 10^(-3)
writedata_steps = 100
T = 100
sigma_to_R = 1

enable_pbc = true

N = N1*N2
time_steps = trunc(Int,T/dt)
M = trunc(Int,T/(dt*writedata_steps)) + 1

for Pe in [0,20]
    traj = run_simulation(Pe, enable_pbc=enable_pbc)
    plot_trajectories(traj, Pe)
end