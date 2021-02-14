module system_methods

include("structs.jl")
using LinearAlgebra, DelimitedFiles

export new_system, evolve_system!, simulate_until_equilibrium!, average_kinetic_energy, average_speed, simulate_and_log_avg_kinetic


const VERTICAL_WALL = -1
const HORISONTAL_WALL = -2


function new_system(
        # Construct new System
        x_positions::Array{Float64},
        y_positions::Array{Float64},
        x_velocities::Array{Float64},
        y_velocities::Array{Float64},
        masses::Array{Float64},
        radii::Array{Float64}
        ;
        restitution_coeff::Float64=1.,
        x_length::Float64=1.,
        y_length::Float64=1.,
        dt_log=1.,
    )
    # Construct new, fully initialized system
    num_particles = length(x_positions)
    sys = System(
        num_particles,
        x_positions,
        y_positions,
        x_velocities,
        y_velocities,
        masses,
        radii,
        x_length,
        y_length,
        restitution_coeff,
        dt_log,

        zeros(num_particles),
        PriorityQueue{Collision, Float64}(),
        0.,
        0,
        0,
    )
    # Make sure all particles are within bounds (not overlap check to come maybe)
    for i in 1:num_particles
        if sys.x_positions[i] - sys.radii[i] < 0.
            raise(error("Particle out of bounds -x"))
        elseif sys.y_positions[i] - sys.radii[i] < 0.
            raise(error("Particle out of bounds -y"))
        elseif sys.x_positions[i] + sys.radii[i] > sys.x_length
            raise(error("Particle out of bounds +x"))
        elseif sys.y_positions[i] + sys.radii[i] > sys.y_length
            raise(error("Particle out of bounds +y"))
        end
    end
    find_and_enqueue_all_col!(sys)
    return sys
end


function find_vertical_wall_col(sys::System, i::Int64)
    # Compute if/when a particle collides with a vertical wall
    x = sys.x_positions[i]
    v_x = sys.x_velocities[i]
    r = sys.radii[i]
    if v_x > 0
        return (sys.x_length - r - x) / v_x
    elseif v_x < 0
        return (r - x) / v_x
    else
        return Inf
    end
end


function find_horisontal_wall_col(sys::System, i::Int64)
    # Compute if/when a particle collides with a horisontal wall
    y = sys.y_positions[i]
    v_y = sys.y_velocities[i]
    r = sys.radii[i]
    if v_y > 0
        return (sys.y_length - r - y) / v_y
    elseif v_y < 0
        return (r - y) / v_y
    else
        return abs(Inf)
    end
end


function find_particle_particle_collision(sys::System, i::Int64, j::Int64)
    # Compute if/when two particles collides with a each other
    ri = sys.radii[i]
    rj = sys.radii[j]
    xi = sys.x_positions[i]
    yi = sys.y_positions[i]
    xj = sys.x_positions[j]
    yj = sys.y_positions[j]
    vxi = sys.x_velocities[i]
    vyi = sys.y_velocities[i]
    vxj = sys.x_velocities[j]
    vyj = sys.y_velocities[j]
    dx = [xj - xi, yj - yi]
    dv = [vxj - vxi, vyj - vyi]
    d = (dot(dv, dx))^2 - dot(dv, dv) * (dot(dx, dx) - (ri+rj)^2)
    if (dot(dv, dx) >= 0) | (d <= 0)
        return Inf
    else
        dt = -(dot(dv, dx)+sqrt(d)) / dot(dv, dv)
        if dt < 0
            #=
            This sometimes occurs, seemingly between very close particles.
            println("neg. dt part.: ", dt)
            println(xi, " ", yi)
            println(xj, " ", yj)
            println(vxi, " ", vyi)
            println(vxj, " ", vyj)
            raise(error("neg. particle time"))
            =#
            return Inf
        end
        return dt
    end
end


function vertical_wall_col!(sys::System, i::Int64)
    # Update velocity and collision count of particle after vertical wall collision
    sys.x_velocities[i] *= -sys.restitution_coeff
    sys.y_velocities[i] *= sys.restitution_coeff
    sys.collision_counts[i] += 1
    return nothing
end


function horisontal_wall_col!(sys::System, i::Int64)
    # Update velocity and collision count of particle after vertical wall collision
    sys.x_velocities[i] *= sys.restitution_coeff
    sys.y_velocities[i] *= -sys.restitution_coeff
    sys.collision_counts[i] += 1
    return nothing
end


function particle_particle_col!(sys::System, i::Int64, j::Int64)
    # Upate velocities and collision count of two particles involved in a particle-particle collision
    ri = sys.radii[i]
    rj = sys.radii[j]
    mi = sys.masses[i]
    mj = sys.masses[j]
    xi = sys.x_positions[i]
    yi = sys.y_positions[i]
    xj = sys.x_positions[j]
    yj = sys.y_positions[j]
    vxi = sys.x_velocities[i]
    vyi = sys.y_velocities[i]
    vxj = sys.x_velocities[j]
    vyj = sys.y_velocities[j]
    res_coeff = sys.restitution_coeff
    dx = [xj - xi, yj - yi]
    dv = [vxj - vxi, vyj - vyi]
    vi = [vxi, vyi]
    vj = [vxj, vyj]
    vi = vi + ((1+res_coeff) * (mj/(mi+mj)) * (dot(dv, dx)/((ri+rj)^2))) * dx
    vj = vj - ((1+res_coeff) * (mi/(mi+mj)) * (dot(dv, dx)/((ri+rj)^2))) * dx

    sys.x_velocities[i], sys.y_velocities[i] = vi
    sys.x_velocities[j], sys.y_velocities[j] = vj

    sys.collision_counts[i] += 1
    sys.collision_counts[j] += 1
    sys.num_pp_collisions += 1
    return nothing
end


function collision!(sys::System, col::Collision)
    # Update velocities of a system's particle(s) involved in the given collision
    if col.collider2 > 0
        # positive index i.e. particle
        particle_particle_col!(sys, col.collider1, col.collider2)
    elseif col.collider2 == VERTICAL_WALL
        vertical_wall_col!(sys, col.collider1)
    elseif col.collider2 == HORISONTAL_WALL
        horisontal_wall_col!(sys, col.collider1)
    end
    sys.num_collisions += 1
    return nothing
end


function find_and_enqueue_all_col!(sys::System)
    # Compute and store all next collitions for every particle (used for initialization)
    for i in 1:sys.num_particles
        # Find collision time deltas
        dt_vertical_wall = find_vertical_wall_col(sys, i)
        dt_horisontal_wall = find_horisontal_wall_col(sys, i)
        dt_particle = Inf
        j = 0
        for k in i+1:sys.num_particles
            dtp = find_particle_particle_collision(sys, i, k)
            if dtp < dt_particle
                dt_particle = dtp
                j = k
            end
        end

        # enqueue possible collisions
        if dt_vertical_wall != Inf
            # enque vertical wall collision
            t_vertical = sys.current_time + dt_vertical_wall
            enqueue!(
                sys.collision_queue,
                Collision(
                    t_vertical,
                    i,
                    VERTICAL_WALL,
                    sys.collision_counts[i],
                    0,
                    ),
                t_vertical,
            )
        end
        if dt_horisontal_wall != Inf
            # enque horisontal wall collision
            t_horisontal = sys.current_time + dt_horisontal_wall
            enqueue!(
                sys.collision_queue,
                Collision(
                    t_horisontal,
                    i,
                    HORISONTAL_WALL,
                    sys.collision_counts[i],
                    0,
                    ),
                t_horisontal,
            )
        end
        if dt_particle != Inf
            # enque vertical wall collision
            t_particle = sys.current_time + dt_particle
            enqueue!(
                sys.collision_queue,
                Collision(
                    t_particle,
                    i,
                    j,
                    sys.collision_counts[i],
                    sys.collision_counts[j],
                    ),
                t_particle,
            )
        end
    end
    return nothing
end


function find_and_enqueue_col!(sys::System, i::Int64)
    # Compute and upcoming possible collisions for the given particle i

    # Get collision time deltas
    dt_vertical_wall = find_vertical_wall_col(sys, i)
    dt_horisontal_wall = find_horisontal_wall_col(sys, i)
    dt_particle = Inf
    j = 0
    for k in 1:i-1
        dtp = find_particle_particle_collision(sys, i, k)
        if dtp < dt_particle
            dt_particle = dtp
            j = k
        end
    end
    for k in i+1:sys.num_particles
        dtp = find_particle_particle_collision(sys, i, k)
        if dtp < dt_particle
            dt_particle = dtp
            j = k
        end
    end

    # enqueue possible collisions
    if dt_vertical_wall != Inf
        # enque vertical wall collision
        t_vertical = sys.current_time + dt_vertical_wall
        enqueue!(
            sys.collision_queue,
            Collision(
                t_vertical,
                i,
                VERTICAL_WALL,
                sys.collision_counts[i],
                0,
                ),
            t_vertical,
        )
    end
    if dt_horisontal_wall != Inf
        # enque horisontal wall collision
        t_horisontal = sys.current_time + dt_horisontal_wall
        enqueue!(
            sys.collision_queue,
            Collision(
                t_horisontal,
                i,
                HORISONTAL_WALL,
                sys.collision_counts[i],
                0,
                ),
            t_horisontal,
        )
    end
    if dt_particle != Inf
        # enque vertical wall collision
        t_particle = sys.current_time + dt_particle
        enqueue!(
            sys.collision_queue,
            Collision(
                t_particle,
                i,
                j,
                sys.collision_counts[i],
                sys.collision_counts[j],
                ),
            t_particle,
        )
    end
    return nothing
end


function evolve_system!(sys::System)
    # Evolve system until next collision/event
    # First find next valid collision
    valid = false
    next_col = undef
    while !valid
        next_col = dequeue!(sys.collision_queue)
        if sys.collision_counts[next_col.collider1] == next_col.collision_count1
            if next_col.collider2 < 1
                valid = true
            elseif sys.collision_counts[next_col.collider2] == next_col.collision_count2
                valid = true
            end
        end
    end

    # Translate particles until the time of the first collision
    dt = next_col.time - sys.current_time
    sys.x_positions += dt * sys.x_velocities
    sys.y_positions += dt * sys.y_velocities

    # Translate time forward to the time of the collision
    sys.current_time += dt

    # Update velocity of colliding particle(s)
    collision!(sys, next_col)

    # compute and store new next collision for the particle(s) involved in the collision
    find_and_enqueue_col!(sys, next_col.collider1)
    if next_col.collider2 > 0
        find_and_enqueue_col!(sys, next_col.collider2)
    end
    return nothing
end


function equilibrium_reached(sys::System, mfactor::Int64)
    return sys.num_pp_collisions / sys.num_particles > mfactor
end


function simulate_until_equilibrium!(sys::System; mfactor=100)
    while ! equilibrium_reached(sys, mfactor)
        evolve_system!(sys)
    end
    println("num collisions: ", sys.num_collisions)
end


function simulate_until_equilibrium!(sys::System, max_iter::Int64; mfactor=100)
    i = 0
    while ! equilibrium_reached(sys, mfactor)
        evolve_system!(sys)
        i += 1
        if i == max_iter
            println("max iter reached")
            break
        end
    end
    println("num collisions: ", sys.num_collisions)
end


function average_speed(sys::System; start=1, stop=0)
    i, j = start, stop
    if j == 0
        j = sys.num_particles
    end
    return sum((sys.x_velocities[i:j].^2 + sys.y_velocities[i:j].^2)) / (j - i + 1)
end


function average_kinetic_energy(sys::System; start=1, stop=0)
    i, j = start, stop
    if j == 0
        j = sys.num_particles
    end
    return 0.5 * dot((sys.x_velocities[i:j].^2 + sys.y_velocities[i:j].^2), sys.masses[i:j]) / (j - i + 1)
end


function simulate_and_log_avg_kinetic(sys::System, num_interval_col::Int64, fname::String; mfactor=20)
    open(fname*"avg_energies.csv", "w") do output
        writedlm(output, ["m0" "4m0" "sys"], ',')
        while ! equilibrium_reached(sys, mfactor)
            evolve_system!(sys)
            if sys.num_particles % num_interval_col == 0
                avg_kinetic_energy_m0 = average_kinetic_energy(sys, stop=Int(sys.num_particles/2))
                avg_kinetic_energy_4m0 = average_kinetic_energy(sys, start=Int(sys.num_particles/2))
                avg_kinetic_energy_all = average_kinetic_energy(sys)
                writedlm(output, [avg_kinetic_energy_m0 avg_kinetic_energy_4m0 avg_kinetic_energy_all], ',')
            end
        end
    end
end


end
