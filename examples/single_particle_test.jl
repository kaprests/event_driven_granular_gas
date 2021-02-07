println("Getting PyPlot...")
using PyPlot
println("Done getting PyPlot")

include("../src/granular_gas.jl")
using .granulargas


if PROGRAM_FILE == basename(@__FILE__)
    # create simple system
    println("creating system")
    N = 1
    r = 1.
    Lx = 10.
    Ly = 10.
    system = new_system(
        [5.],
        [8.],
        [15.],
        [15.],
        [1.],
        [r],
        x_length=Lx,
        y_length=Ly,
    )

    num_col = 10000
    x = zeros(num_col+1)
    y = zeros(num_col+1)
    x[1] = system.x_positions[1]
    y[1] = system.y_positions[1]
    for i in 1:num_col
        # evolve
        evolve_system!(system)
        x[i+1] = system.x_positions[1]
        y[i+1] = system.y_positions[1]
    end
    # plot path
    xlim(0, Lx)
    ylim(0, Ly)
    plot(x, y)
    show()
end
