println("Getting PyPlot...")
using PyPlot
println("Done getting PyPlot")

include("../../granular_gas.jl")
using .granulargas

const PLOT = true

# system config
const num_particles = 5001 # number of particles

# Wall particles
const r = 0.004 # pack fraction aprox .51 for 5000 particles in the lower half
const m = 100. # mass (common for all)
const Lx = 1. # Length of box in x-dir
const Ly = 1. # Length of box in y-dir
const Nx = 150 # Discretization number x-dir for initial positions
const Ny = 75 # Discretization number y-dir for initial positions
x_pos, y_pos = random_positions(num_particles, r, Lx, 0.5*Ly, Nx, Ny)
x_vel, y_vel = zeros(num_particles), zeros(num_particles)
masses = ones(num_particles)*m
radii = ones(num_particles)*r

# Projectile
const x0 = 0.5
const y0 = 0.75
const v0_x = 0.
const v0_y = -1.
const m_proj = 25*m
const r_proj = 5*r
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
    restitution_coeff=0.5,
)
println("System created!")

# save snapshot data
save_pos(system, "initial_task4")

if PLOT
    # Plot initial pos and velocities of system
    xlim(0, Lx)
    ylim(0, Ly)
    scatter(system.x_positions[1], system.y_positions[1], s=125.)
    scatter(system.x_positions[2:end], system.y_positions[2:end], s=5.)
    show()
end

initial_energy = total_energy(system)
target_energy = initial_energy * 0.1

# simulate
println("simulating ")
while target_energy < total_energy(system)
    evolve_system!(system)
end

println("Done simulating")

# save snaphot
save_pos(system, "final_task4")

if PLOT
    # Plot final pos and velocities of system
    title("Final positions")
    xlim(0, Lx)
    ylim(0, Ly)
    scatter(system.x_positions[1], system.y_positions[1], s=125.)
    scatter(system.x_positions[2:end], system.y_positions[2:end], s=5.)
    show()
end
