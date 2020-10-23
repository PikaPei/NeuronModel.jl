using Plots
include("LIF.jl")

neu = LIF("neu")
potential = simulate(neu, neu.rest, 1.0, 100, 0.1)

p = plot(potential, legend=false, dpi=200)
savefig(p, "neu.png")
