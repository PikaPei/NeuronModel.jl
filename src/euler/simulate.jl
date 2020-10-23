using Plots
include("network.jl")

net = Network()
add_neuron(net, "neu1")
add_neuron(net, "neu2")

add_receptor(net, "AMPA", 2.0, 0.0)
add_receptor(net, "GABA", 5.0, -90.0)
set_neuron_receptor_all(net, "AMPA", "GABA")

add_synapse(net, "neu1", "neu2", "AMPA", 25.0)

potential = simulate(net, 100, 0.1)

p = plot(potential, legend=false, dpi=200)
savefig(p, "neu.png")
