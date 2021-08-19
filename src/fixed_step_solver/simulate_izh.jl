using Plots
include("network.jl")

net = Network()

add_neuron(net, "neux", Izhikevich)
add_neuron(net, "neuy", Izhikevich)

add_receptor(net, "AMPA", 2.0, 0.0)
add_receptor(net, "GABA", 5.0, -90.0)
set_neuron_receptor_all(net, "AMPA", "GABA")

add_synapse(net, "neux", "neuy", "AMPA", 10.0)

add_event(net, 0.1, "Current", "neux", 70.0, 10.0)
add_event(net, 0.1, "Current", "neuy", 0.0, 10.0)
add_event(net, 1000.0, "End")

results = simulate(net, Izhikevich; dt=0.1, store_potential=true)
# results = simulate(net, Izhikevich; solver=RK4(), dt=0.1, store_potential=true)

p = plot(results[:, 2:end], legend=false, dpi=200)
savefig(p, "neu.png")
