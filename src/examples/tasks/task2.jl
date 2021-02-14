println("Getting PyPlot...")
using PyPlot
println("Done getting PyPlot")

include("../../granular_gas.jl")
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

# Histogram of initial velocities
hist(sqrt.(system.x_velocities[1:Int(system.num_particles/2)].^2
    + system.y_velocities[1:Int(system.num_particles/2)].^2), label="m0")
hist(sqrt.(system.x_velocities[Int(system.num_particles/2):end].^2
    + system.y_velocities[Int(system.num_particles/2):end].^2), label="4m0")
show()

num_snapshots = 5
println("simulating, number of equilibrium snapshots: ", num_snapshots)
for i in 1:num_snapshots
    simulate_until_equilibrium!(system, system.num_particles*100, mfactor=10)
    hist(sqrt.(system.x_velocities[1:Int(system.num_particles/2)].^2
        + system.y_velocities[1:Int(system.num_particles/2)].^2), label="mass: m0"*"snapshot #: "*string(i))
    hist(sqrt.(system.x_velocities[Int(system.num_particles/2):end].^2
        + system.y_velocities[Int(system.num_particles/2):end].^2), label="mass: 4m0"*"snapshot #: "*string(i))
    system.num_pp_collisions = 0
    system.num_collisions = 0
end
println("Done simulating")
legend()
show()

# Average speed and kinetic energy for particles
avg_speed_m0 = average_speed(system, stop=Int(system.num_particles/2))
avg_speed_4m0 = average_speed(system, start=Int(system.num_particles/2))
avg_kinetic_energy_m0 = average_kinetic_energy(system, stop=Int(system.num_particles/2))
avg_kinetic_energy_4m0 = average_kinetic_energy(system, start=Int(system.num_particles/2))
println("average speed m0: ", avg_speed_m0, " | average speed 4m0: ", avg_speed_4m0)
println("average kinetic energy m0: ", avg_kinetic_energy_m0, " | average kinetic energy 4m0: ", avg_kinetic_energy_4m0)

# Plot final pos and velocities of system
title("Final positions")
scatter(system.x_positions, system.y_positions, s=5.)
show()
title("Final velocities")
scatter(system.x_velocities, system.y_velocities, s=5.)
show()
