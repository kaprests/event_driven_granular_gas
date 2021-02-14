module utils

using StatsBase
using PyPlot

export random_positions, random_directions


function random_positions(N::Int64, r::Float64, Lx::Float64, Ly::Float64, Nx::Int64, Ny::Int64)
    # Generates ish random positions for N particles
    r *= 1.5 # avoid particles starting in collision
    if ((Lx-2*r)/(Nx-1) < r) | ((Ly-2*r)/(Ny-1) < r)
        return error("Invalid parameters, cannot guarantee no overlap")
    elseif N > Nx*Ny
        return error("Invalid parameters, cannot draw N positions")
    end
    xy_grid = CartesianIndices(zeros(Nx, Ny))
    rand_grid_idx = sample(1:Nx*Ny, N, replace=false)
    xy_idx = xy_grid[rand_grid_idx]
    x_idx = getindex.(xy_idx, 1)
    y_idx = getindex.(xy_idx, 2)
    x_pos = r .+ (x_idx .- 1) * (Lx-2*r)/(Nx-1)
    y_pos = r .+ (y_idx .- 1) * (Ly-2*r)/(Ny-1)
    return x_pos, y_pos
end


function random_directions(N::Int64)
    # generate random directions
    angles = rand(0:0.01:2*pi, N)
    norm_complex = exp.(angles*im)
    return real(norm_complex), imag(norm_complex) # vx, vy
end


if PROGRAM_FILE == basename(@__FILE__)
    # run (eventual) tests
    println("done importing, running function")

    # test random_positions
    x, y = random_positions(10000, 1e-9, 1., 1., 1000, 1000)
    title("random positions")
    scatter(x, y, s=2)
    show()

    # test random_directions
    title("random directions")
    x, y = random_directions(1000)
    scatter(x, y, s=2)
    show()
end


end
