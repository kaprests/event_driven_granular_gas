module granulargas

include("system_methods.jl")
include("utils.jl")

using .system_methods, .utils
export random_positions, random_directions, new_system, evolve_system!
end
