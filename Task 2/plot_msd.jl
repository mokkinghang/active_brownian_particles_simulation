include("Task2.jl")

N1, N2 = 40,25
N = N1*N2
L = 85
dt = 10^(-3)
writedata_steps = 10
T = 100
L = 85 # side length of box
dt = 10^(-3) # time step
writedata_steps = 10 # number of time steps to write data into matrix
mu = 1 # mobility
h = 1 # hardness of sphere
R = 2^(1/6)  # Cut off radius

time_steps = trunc(Int,T/dt)
M = trunc(Int,T/(dt*writedata_steps)) + 1
sigma_to_R = 1/R # The ratio sigma/R
number_of_cells = trunc(Int,L/R)  # Number of cells ACROSS ONE DIMENSION, the actual no. of cells in the entire domain is the square of it
cell_size = L/number_of_cells

enable_pbc = true
min_img_convention = true
cell_linked_list = true

Pe_arr = [0,5,10,15,20]

save_traj_to_jld2(Pe_arr; enable_pbc=enable_pbc, min_img_convention=min_img_convention, cell_linked_list=cell_linked_list)
plot_msd(["Pe=0.jld2", "Pe=5.jld2", "Pe=10.jld2", "Pe=15.jld2", "Pe=20.jld2"])
