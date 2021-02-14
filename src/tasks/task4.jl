println("Getting PyPlot...")
using PyPlot
println("Done getting PyPlot")

include("../granular_gas.jl")
using .granulargas


# system config
num_particles = 5001 # number of particles

# Wall particles
r = 0.0025 # radius of (same for every)
m = 1 # mass (common for all)
Lx = 1. # Length of box in x-dir
Ly = 1. # Length of box in y-dir
Nx = 200 # Discretization number x-dir for initial positions
Ny = 100 # Discretization number y-dir for initial positions
x_pos, y_pos = random_positions(num_particles, r, Lx, 0.5*Ly, Nx, Ny)
x_vel, y_vel = zeros(num_particles), zeros(num_particles)
masses = ones(num_particles)*m
radii = ones(num_particles)*r

# Projectile
x0 = 0.5
y0 = 0.75
v0_x = 0.
v0_y = -1.
m_proj = 5.
r_proj = 0.
x_pos[1] = x0
y_pos[1] = y0
x_vel[1] = v0_x
y_vel[1] = v0_y
masses[1] = m_proj
radii[1] = r_proj


println("creating system")
system = new_system(
    x_pos,
    y_pos,
    x_vel,
    y_vel,
    masses,
    radii,
    x_length=Lx,
    y_length=Ly,
)
println("System created!")

# Plot initial pos and velocities of system
xlim(0, Lx)
ylim(0, Ly)
scatter(system.x_positions, system.y_positions, s=5.)
show()

# simulate
N_col = 5000
println("simulating ", N_col, " collisions")
for i in 1:N_col
    evolve_system!(system)
end
println("Done simulating")

# Plot final pos and velocities of system
title("Final positions")
xlim(0, Lx)
ylim(0, Ly)
scatter(system.x_positions, system.y_positions, s=5.)
show()
