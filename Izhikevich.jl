using DataStructures


abstract type NeuralPopulation end


mutable struct Izhikevich <: NeuralPopulation
	name::String
	N::Int
	C::Float64
	Vr::Float64
	Vth::Float64
	Vpeak::Float64
	k::Float64
	a::Float64
	b::Float64
	c::Float64
	d::Float64
end

function Izhikevich(name;
					N=1, C=100.0, Vr=-60.0, Vth=-40.0, Vpeak=35.0,
					k=0.7, a=0.03, b=-2.0, c=-50.0, d=100.0)
	return Izhikevich(name, N, C, Vr, Vth, Vpeak, k, a, b, c, d)
end


mutable struct Network
	neu::OrderedDict{String, NeuralPopulation}

	Network() = new(OrderedDict{String,NeuralPopulation}())
end

# function Base.show(io::IO, net::Network)
#     println(io, keys(net.neu))
# end

default_model = Izhikevich
function add_neuron(net::Network, name, model=default_model; kwargs...)
	neuron = default_model(name; kwargs...)
	net.neu[name] = neuron
end
