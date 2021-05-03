println("Getting PyPlot...")
using PyPlot
println("Done getting PyPlot")
include("../../granular_gas.jl")
using .granulargas


## system configuration ###

num_particles = 2000 # number of particles
r = 0.0025 # radius of (same for every)
v0 = .1 # initial velocity (magnitude)
m = 1e-1 # mass (common for all)
Lx = 1. # Length of box in x-dir
Ly = 1. # Length of box in y-dir
Nx = 100 # Discretization number x-dir for initial positions
Ny = 100 # Discretization number y-dir for initial positions

# Initial positions and velocities
x_pos, y_pos = random_positions(num_particles, r, Lx, Ly, Nx, Ny)
x_vel, y_vel = random_directions(num_particles)
x_vel *= v0
y_vel *= v0

## Create system ##

println("creating system")
system = new_system(
    # required params
    x_pos,
    y_pos,
    x_vel,
    y_vel,
    ones(num_particles)*m,
    ones(num_particles)*r,
    
    # optional params
    x_length=Lx,
    y_length=Ly,
)
println("System created!")

# Histogram of initial velocities
hist(sqrt.(system.x_velocities.^2 + system.y_velocities.^2))
show()

## simulation ##
num_snapshots = 5
println("simulating, number of equilibrium snapshots: ", num_snapshots)
for i in 1:num_snapshots
    simulate_until_equilibrium!(system, system.num_particles*100, mfactor=10)
    hist(sqrt.(system.x_velocities.^2 + system.y_velocities.^2), label=string(i))
    system.num_pp_collisions = 0
    system.num_collisions = 0
end
println("Done simulating")
legend()
show()

# Plot final pos and velocities of system
title("Final positions")
scatter(system.x_positions, system.y_positions, s=5.)
show()
title("Final velocities")
scatter(system.x_velocities, system.y_velocities, s=5.)
show()
