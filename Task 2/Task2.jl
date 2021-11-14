"""
Parameters description:
N = N1*N2: number of particles, where N1 is number of particles along x-direction and N2 along y-direction
Pe: Peclet number
L: side length of box
dt: time step size
writedata_steps: number of time steps to write data into matrix
T: Total simulation time
sigma_to_R: The ratio sigma/R, where R is the radius of particle (cut-off radius at the same time)
time_steps: Number of time steps
M: number of "screenshots" made, i.e. writing data to matrix
mu: particle mobility
h: hardness of particle
R: radius particle, also the cut-off radius (often taken as 2^(1/6)) which is the length of cell in cell linked list implementation
number_of_cells: Number of cell in one dimension (Total number of cells in the entire domain is its square)
cell_size: side length of square cell, which is not always the same as R, since we have to make number_of_cells an integer. It is mostly slightly larger than R.
enable_pbc: boolean value, whether periodic boundary condition is enabled
min_img_convention: boolean value, whether minimum image convention is enabled
cell_linked_list: boolean value, whether cell linked list is implemented for pair interactions
"""

using Random, Plots, LinearAlgebra, Statistics, LaTeXStrings, JLD2
pyplot()

function update!(arr, Pe; enable_pbc=true, min_img_convention=true, cell_linked_list=true)
    # The displacement in one time step due to self-propulsion and transform it to 3xN matrix
    self_prop_term = [Pe*cos.(arr[3,:]), Pe*sin.(arr[3,:]), zeros(N)] * dt
    self_prop_disp = transpose(hcat(self_prop_term...))
    # The displacement in one time step due to stochasticity
    stochastic_term = [sqrt(2), sqrt(2), sqrt(3/2)*sigma_to_R] * sqrt(dt)
    stochastic_disp = repeat(stochastic_term, outer = [1,N]) .* randn(3,N)
    # Get pair forces
    forces = get_pair_forces(arr; enable_pbc=enable_pbc, min_img_convention=min_img_convention, cell_linked_list=cell_linked_list)
    # The displacement due to particle pair interactions
    force_disp = forces * mu * dt
    # Update the coordinates
    arr[:,:] += self_prop_disp + stochastic_disp + force_disp
    # Periodic boundary condition 
    if enable_pbc
        pbc!(arr)
    end
    return arr
end

function run_simulation(Pe; enable_pbc=true, min_img_convention=true, cell_linked_list=true)
    """
    We first define a 3*N*M multidimensional array, representing 3 degrees of freedom, N particles and M screen shots.
    To save memory, not every time step will be saved into the matrix. 
    """
    traj = zeros(3,N,M)
    step_count = 0
    entry_count = 1

    arr = even_spacing_initialization()
    traj[:,:,1] = arr

    for t in 1:time_steps
        update!(arr, Pe; enable_pbc=enable_pbc, min_img_convention=min_img_convention, cell_linked_list=cell_linked_list)
        step_count += 1
        if step_count % writedata_steps == 0
            entry_count += 1
            traj[:,:,entry_count] = arr
        end
    end
    return traj
end

function save_traj_to_jld2(Pe_arr; enable_pbc=true, min_img_convention=min_img_convention, cell_linked_list=cell_linked_list)
    for Pe in Pe_arr
        println("Running Pe=$(Pe)...")
        traj = run_simulation(Pe; enable_pbc=enable_pbc, min_img_convention=min_img_convention, cell_linked_list=cell_linked_list)
        t_arr = [dt*writedata_steps*(t-1) for t in 1:M]
    
        jldopen("Pe=$(Pe).jld2", "w") do file
            file["Pe"] = Pe
            file["traj"] = traj
            file["t_arr"] = t_arr
            file["msd"] = msd(traj)
        end
    end
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
    x_arr = vcat([[p*i+p/2 for i in 0:N1-1] for i in 0:N2-1]...)
    y_arr = vcat([[q*j+q/2 for i in 0:N1-1] for j in 0:N2-1]...)
    arr[1,:] = x_arr
    arr[2,:] = y_arr
    arr[3,:] = rand(1,N) * 2*pi
    return arr
end

function draw_particles(arr; markersize=10)
    x_arr = arr[1,:]
    y_arr = arr[2,:]
    colors_list = ["red","blue","yellow","orange","brown","green","purple","black","gray"]
    colors = vcat(vcat([colors_list for i in 1:trunc(Int,N/size(colors_list)[1])]...),
    colors_list[1:trunc(Int,N%size(colors_list)[1])])
    scatter(x_arr, y_arr, legend=false, xlims=(0,L), ylims=(0,L), markersize=markersize, markercolor=colors)
end

function mirror_images(pos)
    """
    Given position as a 2-element vector, return an array of positions of all mirror images
    """
    x,y = Tuple(pos)
    return [[x,y],[x+L,y],[x-L,y],[x,y+L],[x,y-L],[x+L,y+L],[x+L,y-L],[x-L,y+L],[x-L,y-L]]
end

function get_dist(pos_i, pos_j)
    """
    A function returning the distance, x-y components of joining lines of two particles
    Input: two 2-element vectors
    """
    j_mirror_images = mirror_images(pos_j)
    dist_arr = [norm(pos_i-img_pos_j) for img_pos_j in j_mirror_images]
    min_dist, min_dist_mirror_image_j = findmin(dist_arr)[1], j_mirror_images[findmin(dist_arr)[2]]
    vector_ij = min_dist_mirror_image_j - pos_i
    dx, dy = vector_ij[1], vector_ij[2]
    return min_dist, dx, dy
end

function get_dist_2(pos_i, pos_j)
    """
    When min_img_convention is disabled
    """
    dist = norm(pos_i-pos_j)
    vector_ij = pos_j - pos_i
    dx, dy = vector_ij[1], vector_ij[2]
    return dist, dx, dy
end

# Potential force (L-J potential)
F(r) = h*(48*r^(-13) - 24*r^(-7))

function get_pair_forces(arr; enable_pbc=true, min_img_convention=true, cell_linked_list=true)
    """
    Return forces as 3 x N matrix. The third row will be zeros if there is no (deterministic) external torque.
    """
    forces = zeros(3,N)
    interaction_pairs = cell_linked_list ? get_interaction_pairs(arr) : [(i,j) for i in 1:N, j in 1:N if i<j]
    for (i,j) in interaction_pairs
        pos_i, pos_j = arr[1:2,i], arr[1:2,j]
        r,dx,dy = min_img_convention ? get_dist(pos_i, pos_j) : get_dist_2(pos_i, pos_j)
        force = F(r)
        force_vec = [force*dx/r, force*dy/r, 0]
        forces[:,i] += -force_vec
        forces[:,j] += force_vec
    end
    return forces
end

function new_cell_dict()
    """
    Create a dictionary whose keys are all possible cells, while the values are empty arrays
    The values will later be the arrays containing particles situating the cells corresponding to the keys
    """
    tuples_arr = [[(i,j) for i in 1:number_of_cells, j in 1:number_of_cells]...]
    empty_cell_dict = Dict([(tuple,[]) for tuple in tuples_arr])
    return empty_cell_dict
end

function get_cell_dict_and_array(arr)
    """
    Return a dictionary cell_dict. Keys: cell; Values: particles in the cell.
    Return an array cell_arr. cell_arr[i] gives the cell at which particle i is located in.
    """
    cell_dict = new_cell_dict()  # Create a dictionary 
    x_arr, y_arr = arr[1,:], arr[2,:]
    cell_x,cell_y = trunc.(Int,x_arr./cell_size).+1,trunc.(Int,y_arr./cell_size).+1
    cell_arr = [i for i in zip(cell_x,cell_y)]
    for (i,cell) in enumerate(cell_arr)
        push!(cell_dict[cell],i) 
    end
    return cell_dict, cell_arr
end

function get_neighboring_cells(cell)
    """
    Input: cell (Tuple{Int64, Int64})
    Return: array containing tuples representing the neighboring cells and the original cell
    """
    @assert typeof(cell) == Tuple{Int64, Int64}
    @assert cell[1] > 0 && cell[1] <= number_of_cells
    @assert cell[2] > 0 && cell[2] <= number_of_cells
    cell_x,cell_y = cell
    neighboring_cells_raw = [
    (cell_x-1,cell_y-1),(cell_x-1,cell_y),(cell_x-1,cell_y+1),
    (cell_x,cell_y-1),(cell_x,cell_y),(cell_x,cell_y+1),
    (cell_x+1,cell_y-1),(cell_x+1,cell_y),(cell_x+1,cell_y+1)
    ]
    neighboring_cells = [mod.(raw_cell.-1,number_of_cells).+1 for raw_cell in neighboring_cells_raw]
    return neighboring_cells
end

function get_interaction_pairs(arr)
    """
    Return an array that contains all tuples indicating interaction pairs
    e.g. (4,7) indicates that particles 4 and 7 will interact
    No interaction if the pair not in array
    """
    interaction_pairs = []
    cell_dict, cell_arr = get_cell_dict_and_array(arr)
    for i in 1:N
        c1 = cell_arr[i]  # Cell tuple for particle i
        neighboring_cells = get_neighboring_cells(c1) # Arr of all neighboring cells of cell_i
        particle_j_arr = []  # An array that will contain indices of particles interacting with i 
        for c2 in neighboring_cells
            particle_j_arr = vcat(particle_j_arr, [j for j in cell_dict[c2] if j > i])
        end
        interaction_pairs_with_i = [(i,j) for j in particle_j_arr]
        interaction_pairs = vcat(interaction_pairs, interaction_pairs_with_i)
    end
    return interaction_pairs
end

"""
Theoretical MSD (nondimensionalized) for non-interactive particles
"""
theoretical_msd(t, Pe) = (4 + 8/3*Pe^2*(sigma_to_R)^(-2))*t + 32/9*Pe^2*(sigma_to_R)^(-4)*(exp(-3/4*sigma_to_R^2*t)-1)

function msd(traj)
    init_pos_arr = traj[1:2,:,1]
    msd_arr = zeros(M)
    for t in 1:M
        difference = traj[1:2,:,t] - init_pos_arr
        msd_arr[t] = mean(map(norm, eachslice(difference, dims=2)).^2)
    end
    return msd_arr
end

function plot_msd(jld2_file_arr)
    # Open JLD2 file
    msd_plot = plot(xlabel=L"t", ylabel=L"MSD", legend=:bottomright, 
    dpi=300, xaxis=:log, yaxis=:log, xlims=(dt*writedata_steps, T), ylims=(dt,Inf))

    for jld2_file in jld2_file_arr
        file = jldopen(jld2_file, "r")
        Pe = file["Pe"]
        t_arr = file["t_arr"]
        msd_arr = file["msd"]

        theoretical_msd_arr = [theoretical_msd(t, Pe) for t in t_arr]

        plot!(t_arr, msd_arr, label="Pe=$(Pe) (interactive)")
        plot!(t_arr, theoretical_msd_arr, label="Pe=$(Pe) (non-interactive)")
    end

    savefig(msd_plot, "msd.png")
end