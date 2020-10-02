using DataStructures

include("Izhikevich.jl")


"""
Event
"""
abstract type AbstractEvent end

mutable struct EventCurrentInjection <: AbstractEvent
    time::Float64
    type::String
    population::String
    Gauss_mean::Float64
    Gauss_std::Float64
end

mutable struct EventExtFreq <: AbstractEvent
    time::Float64
    type::String
    population::String
    receptor::String
    freq::Float64
end

mutable struct EventEndTrial <: AbstractEvent
    time::Float64
    type::String
end


"""
Network
"""
mutable struct Network
    neu::OrderedDict{String,NeuronModel}
	receptors::OrderedDict{String,Receptor}
    event::Vector{AbstractEvent}
    event_temp::Vector{EventExtFreq}
    event_end::EventEndTrial

    Network() = new(OrderedDict{String,NeuronModel}(),
					OrderedDict{String,Receptor}(),
                    Vector{AbstractEvent}(),
                    Vector{EventExtFreq}(),
                    EventEndTrial(0.0, "EndTrial"))
end

# function Base.show(io::IO, net::Network)
#     println(io, keys(net.neu))
# end

default_model = Izhikevich
function add_neuron(net::Network, name, model=default_model; kwargs...)
    neuron = default_model(name; kwargs...)
    net.neu[name] = neuron
end


function add_receptor(net::Network,
                      name::String,
                      tau::Real,
                      reversal::Real,
                      external_freq::Real,
                      mean_external_efficacy::Real,
                      mean_external_connection::Real)

    rec = Receptor(name, tau, reversal, external_freq, mean_external_efficacy, mean_external_connection)
    net.receptors[name] = rec
end


function set_neuron_receptor_all(net::Network, args...)
    for neu in net.neu, rec in args
        push!(neu.second.receptors, net.receptors[rec])
    end
end


function add_target(net::Network,
                    pre_syn::String,
                    post_syn::String,
                    target_receptor::String,
                    mean_efficacy::Float64,
                    weight::Float64)

    post_syn in keys(net.neu) || begin @error "No target neural population can be found."; exit(1) end
    net.neu[pre_syn].targets[post_syn] = Target(post_syn, target_receptor, mean_efficacy, weight)
end


function add_event(net::Network, time, type::String, args...)
    if type == "Current"
        event = EventCurrentInjection(time, type, args...)
        push!(net.event, event)

    elseif type == "ExtFreq"
        event = EventExtFreq(time, type, args...)
        push!(net.event_temp, event)

    elseif type == "EndTrial"
        event = EventEndTrial(time, type, args...)
        net.event_end = event
    end
end
