include("Izhikevich.jl")


net = Network()
add_neuron(net, "neu"; N=100)
# add_neuron(net, "neu", Izhikevich; N=100)
