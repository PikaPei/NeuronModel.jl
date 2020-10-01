include("network.jl")
include("simulator.jl")


net = Network()
add_neuron(net, "neu")
# add_neuron(net, "neu"; N=100)
# add_neuron(net, "neu", Izhikevich; N=100)


add_event(net, 0.1, "Current", "neu", 70.0, 0.0)
add_event(net, 1000.0, "EndTrial")


prob, tstops, save_idxs = gen_problem(net)
sol = solve(prob, tstops=tstops, saveat=0.1, save_idxs=save_idxs, abstol=1e-9, reltol=1e-9)


output_mempot("MemPotALL.dat", net, sol)
output_spike("SpikeALL.dat", net)
