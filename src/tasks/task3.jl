println("Getting PyPlot...")
using PyPlot
println("Done getting PyPlot")

include("../granular_gas.jl")
using .granulargas


# system config
num_particles = 200 # number of particles
r = 0.0025 # radius of (same for every)
v0 = .1 # initial velocity (magnitude)
m0 = 1 # reference mass
Lx = 1. # Length of box in x-dir
Ly = 1. # Length of box in y-dir
Nx = 100 # Discretization number x-dir for initial positions
Ny = 100 # Discretization number y-dir for initial positions

x_pos, y_pos = random_positions(num_particles, r, Lx, Ly, Nx, Ny)
x_vel, y_vel = random_directions(num_particles)
x_vel *= v0
y_vel *= v0

m1 = ones(Int(num_particles/2))*m0
m2 = ones(Int(num_particles/2))*4*m0
masses = vcat(m1, m2)

println("creating system")
system = new_system(
    x_pos,
    y_pos,
    x_vel,
    y_vel,
    masses,
    ones(num_particles)*r,
    x_length=Lx,
    y_length=Ly,
)
println("System created!")

simulate_and_log_avg_kinetic(system, Int(round(num_particles/100)), "task3")

# Plot final pos and velocities of system
title("Final positions")
scatter(system.x_positions, system.y_positions, s=5.)
show()
title("Final velocities")
scatter(system.x_velocities, system.y_velocities, s=5.)
show()
