println("Getting PyPlot...")
using PyPlot
println("Done getting PyPlot")

include("../granular_gas.jl")
using .granulargas


if PROGRAM_FILE == basename(@__FILE__)
    # create simple system
    println("creating system")
    r = 0.0001
    m = 0.1
    system = new_system(
        [.25, .75],
        [.75, .25],
        [3., 1.],
        [1., 3.],
        [m, m],
        [r, r],
    )

    num_col = 1000
    x1 = zeros(num_col+1)
    y1 = zeros(num_col+1)
    x2 = zeros(num_col+1)
    y2 = zeros(num_col+1)
    x1[1] = system.x_positions[1]
    y1[1] = system.y_positions[1]
    x2[1] = system.x_positions[2]
    y2[1] = system.y_positions[2]
    for i in 1:num_col
        # evolve
        evolve_system!(system)
        x1[i+1] = system.x_positions[1]
        y1[i+1] = system.y_positions[1]
        x2[i+1] = system.x_positions[2]
        y2[i+1] = system.y_positions[2]
    end
    # plot path
    xlim(0, 1)
    ylim(0, 1)
    plot(x1, y1)
    plot(x2, y2)
    show()
end
