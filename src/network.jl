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
    event::Vector{AbstractEvent}
    event_temp::Vector{EventExtFreq}
    event_end::EventEndTrial

    Network() = new(OrderedDict{String,NeuronModel}(),
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
