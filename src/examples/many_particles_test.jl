println("Getting PyPlot...")
using PyPlot
println("Done getting PyPlot")

include("../granular_gas.jl")
using .granulargas


if PROGRAM_FILE == basename(@__FILE__)
    println("creating system")

    num_particles = 1000
    r = 1e-3
    v0 = .1
    m = 1e-1
    Lx = 1.
    Ly = 1.
    Nx = 100
    Ny = 100

    x_pos, y_pos = random_positions(num_particles, r, Lx, Ly, Nx, Ny)
    x_vel, y_vel = random_directions(num_particles)
    x_vel *= v0
    y_vel *= v0

    system = new_system(
        x_pos,
        y_pos,
        x_vel,
        y_vel,
        ones(num_particles)*m,
        ones(num_particles)*r,
        x_length=Lx,
        y_length=Ly,
    )

    # Plot initial pos and velocities of system
    xlim(0, Lx)
    ylim(0, Ly)
    scatter(system.x_positions, system.y_positions, s=5.)
    show()
    xlim(-1.5*v0, 1.5*v0)
    ylim(-1.5*v0, 1.5*v0)
    scatter(system.x_velocities, system.y_velocities, s=5.)
    show()

    loops = 100*num_particles
    x1 = zeros(loops+1)
    y1 = zeros(loops+1)
    x1[1] = system.x_positions[1]
    y1[1] = system.y_positions[1]
    for i in 1:loops
        evolve_system!(system)
        x1[i+1] = system.x_positions[1]
        y1[i+1] = system.y_positions[1]
    end

    # plot path
    plot(x1, y1)
    plot(x1[1], y1[1], "bo")
    plot(x1[end], y1[end], "^")
    println("start: ", x1[1], " ",  y1[1])
    println("slutt: ", x1[end], " ",  y1[end])
    show()

    # Plot final pos and velocities of system
    title("Final positions")
    scatter(system.x_positions, system.y_positions, s=5.)
    show()
    title("Final velocities")
    scatter(system.x_velocities, system.y_velocities, s=5.)
    show()
end
