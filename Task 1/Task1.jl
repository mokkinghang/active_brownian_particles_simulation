"""
Parameters description:
N = N1*N2: number of particles, where N1 is number of particles along x-direction and N2 along y-direction
Pe: Peclet number
L: side length of box
dt: time step size
writedata_steps: number of time steps to write data into matrix
T: Total simulation time
sigma_to_R: The ratio sigma/R, where R is the size of particle
time_steps: Number of time steps
M: number of "screenshots" made, i.e. writing data to matrix
"""

using Random, Plots, LinearAlgebra, Statistics, LaTeXStrings, JLD2
pyplot()

"""
Theoretical MSD (nondimensionalized)
"""
theoretical_msd(t, Pe) = (4 + 8/3*Pe^2*(sigma_to_R)^(-2))*t + 32/9*Pe^2*(sigma_to_R)^(-4)*(exp(-3/4*sigma_to_R^2*t)-1)

function update!(arr, Pe; enable_pbc=true)
    # Loop through each particle
    for i in 1:N
        arr[:,i] += [Pe*cos(arr[3,i]), Pe*sin(arr[3,i]), 0] * dt + [sqrt(2), sqrt(2), sqrt(3/2)*sigma_to_R] .* randn(3) * sqrt(dt)        
    end
    # Periodic boundary condition
    if enable_pbc
        pbc!(arr)
    end
    return arr
end

function pbc!(arr)
    arr[1:2, :] = mod.(arr[1:2, :], L)
    arr[3, :] = mod.(arr[3, :], 2*pi)
    return arr
end

function random_initialization()
    arr = zeros(3,N)
    arr[1:2,:] = rand(2,N) * L
    arr[3,:] = rand(1,N) * 2*pi
    return arr
end

function even_spacing_initialization()
    arr = zeros(3,N)
    p,q = L/N1, L/N2
    x_arr = vcat([[p*i for i in 0:N1-1] for i in 0:N2-1]...)
    y_arr = vcat([[q*j for i in 0:N1-1] for j in 0:N2-1]...)
    arr[1,:] = x_arr
    arr[2,:] = y_arr
    arr[3,:] = rand(1,N) * 2*pi
    return arr
end

function msd(traj)
    init_pos_arr = traj[1:2,:,1]
    msd_arr = zeros(M)
    for t in 1:M
        difference = traj[1:2,:,t] - init_pos_arr
        msd_arr[t] = mean(map(norm, eachslice(difference, dims=2)).^2)
    end
    return msd_arr
end

function msad(traj)
    init_angle_arr = traj[3,:,1]
    msad_arr = zeros(M)
    for t in 1:M
        difference = traj[3,:,t] - init_angle_arr
        msad_arr[t] = mean(difference.^2)
    end
    return msad_arr    
end

function run_simulation(Pe; enable_pbc=true)
    """
    We first define a 3*N*M multidimensional array, representing 3 degrees of freedom, N particles and M screen shots.
    To save memory, not every time step will be saved into the matrix. 
    """
    traj = zeros(3,N,M)
    step_count = 0
    entry_count = 1
    arr = random_initialization()
    traj[:,:,1] = arr

    for t in 1:time_steps
        update!(arr, Pe, enable_pbc=enable_pbc)
        step_count += 1
        if step_count % writedata_steps == 0
            entry_count += 1
            traj[:,:,entry_count] = arr
        end
    end
    return traj
end

function save_traj_to_jld2(Pe_arr; enable_pbc=true)
    for Pe in Pe_arr
        println("Running Pe=$(Pe)...")
        traj = run_simulation(Pe, enable_pbc=enable_pbc)
        t_arr = [dt*writedata_steps*(t-1) for t in 1:M]
    
        jldopen("Pe=$(Pe).jld2", "w") do file
            file["Pe"] = Pe
            file["traj"] = traj
            file["t_arr"] = t_arr
            file["msd"] = msd(traj)
            file["msad"] = msad(traj)
        end
    end
end

function plot_trajectories(traj, Pe)
    trajectory_plot = plot(xlabel=L"x", ylabel=L"y", dpi=300, title=L"Pe=%$(Pe)")
    shapes = [:circle, :utriangle, :cross, :star5, :pentagon]
    particle_num = size(traj)[2]
    for i in 1:particle_num
        x_arr = traj[1,i,:]
        y_arr = traj[2,i,:]
        scatter!(x_arr, y_arr, markersize=5, markershape=shapes[i], label="particle $(i)")
    end
    savefig(trajectory_plot, "traj_Pe=$(Pe).png")
end

function plot_msd(jld2_file_arr)
    # Open JLD2 file
    msd_plot = plot(xlabel=L"t", ylabel=L"MSD", legend=:bottomright, 
    dpi=300, xaxis=:log, yaxis=:log, xlims=(dt*writedata_steps, Inf), ylims=(dt,Inf))

    for jld2_file in jld2_file_arr
        file = jldopen(jld2_file, "r")
        Pe = file["Pe"]
        t_arr = file["t_arr"]
        msd_arr = file["msd"]

        theoretical_msd_arr = [theoretical_msd(t, Pe) for t in t_arr]

        plot!(t_arr, msd_arr, label="Pe=$(Pe) (simulated)")
        plot!(t_arr, theoretical_msd_arr, label="Pe=$(Pe) (theoretical)")
    end

    savefig(msd_plot, "msd.png")
end

function plot_msad(jld2_file_arr)
    msad_plot = plot(xlabel=L"t", ylabel=L"MSAD", legend=:bottomright, dpi=300)
    
    for jld2_file in jld2_file_arr
        file = jldopen(jld2_file, "r")
        Pe = file["Pe"]
        t_arr = file["t_arr"]
        msad_arr = file["msad"]

        plot!(t_arr, msad_arr, label="Pe=$(Pe)")
    end

    savefig(msad_plot, "msad.png")
end