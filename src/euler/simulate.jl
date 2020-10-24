using Plots
include("network.jl")

net = Network()

add_neuron(net, "neu1")
add_neuron(net, "neu2")

add_receptor(net, "AMPA", 2.0, 0.0)
add_receptor(net, "GABA", 5.0, -90.0)
set_neuron_receptor_all(net, "AMPA", "GABA")

add_synapse(net, "neu1", "neu2", "AMPA", 25.0)

add_event(net, 0.1, "Current", "neu1", 0.68, 0)
add_event(net, 100.0, "End")

results = simulate(net; dt=0.1, store_potential=true)
# results = simulate(net; solver=RK4(), dt=0.1, store_potential=true)

p = plot(results, legend=false, dpi=200)
savefig(p, "neu.png")
