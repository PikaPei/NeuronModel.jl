include("network.jl")
include("simulator.jl")


net = Network()
add_neuron(net, "neux")
add_neuron(net, "neuy")
# add_neuron(net, "neu"; N=100)
# add_neuron(net, "neu", Izhikevich; N=100)


add_receptor(net, "AMPA", 2.0, 0.0, 0.0, 2.1, 1.0)
add_receptor(net, "GABA", 5.0, -90.0, 0.0, 0.0, 0.0)
set_neuron_receptor_all(net, "AMPA", "GABA")


add_target(net, "neux", "neuy", "AMPA", 1.0, 10.0)


add_event(net, 0.1, "Current", "neux", 70.0, 0.0)
add_event(net, 1000.0, "EndTrial")


prob, tstops, save_idxs = gen_problem(net)
sol = solve(prob, tstops=tstops, saveat=0.1, save_idxs=save_idxs, abstol=1e-9, reltol=1e-9)


output_mempot("MemPotALL.dat", net, sol)
output_spike("SpikeALL.dat", net)
