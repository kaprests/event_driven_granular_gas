###############
### Structs ###
###############


using DataStructures


struct Collision
    time::Float64 # time of collision
    collider1::Int64 # particle index
    collider2::Int64 # particle index, may be negative (corresponding to a wall)
    collision_count1::Int64 # from when the collision was computed
    collision_count2::Int64 # from when the collision was computed
end


mutable struct System
    ### user config ###
    num_particles::Int64
    x_positions::Array{Float64}
    y_positions::Array{Float64}
    x_velocities::Array{Float64}
    y_velocities::Array{Float64}
    masses::Array{Float64}
    radii::Array{Float64}
    x_length::Float64
    y_length::Float64
    restitution_coeff::Float64
    dt_log::Float64

    ### Non user config ###
    collision_counts::Array{Int64}
    collision_queue::PriorityQueue{Collision, Float64} # {Collision, collision time}
    current_time::Float64
    num_pp_collisions::Int64
    num_collisions::Int64
end
