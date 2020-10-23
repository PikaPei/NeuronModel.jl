using DataStructures
using Parameters

abstract type AbstractNeuronModel end

include("receptor.jl")
include("synapse.jl")
include("event.jl")
include("LIF.jl")
include("Izhikevich.jl")


mutable struct Network
    neuron::OrderedDict{String,AbstractNeuronModel}
    receptor::OrderedDict{String,Receptor}
    synapse::OrderedDict{String,Vector{Synapse}}
    event::Vector{AbstractEvent}
    event_end::EventEndTrial

    Network() = new(OrderedDict{String,AbstractNeuronModel}(),
                    OrderedDict{String,Receptor}(),
                    OrderedDict{String,Vector{Synapse}}(),
                    Vector{AbstractEvent}(),
                    EventEndTrial(0.0, "EndTrial"))
end

function Base.show(io::IO, net::Network)
    println(io, keys(net.neuron))
end

default_model = LIF
function add_neuron(net::Network, name, model=default_model; kwargs...)
    neuron = default_model(name; kwargs...)
    net.neuron[name] = neuron
end


function add_receptor(net::Network,
                      name::String,
                      tau::Real,
                      reversal::Real)
    rec = Receptor(name, tau, reversal)
    net.receptor[name] = rec
end


function set_neuron_receptor_all(net::Network, args...)
    for neu in values(net.neuron), rec in args
        push!(neu.receptor, net.receptor[rec])
    end
end


function add_synapse(net::Network, pre::String, post::String,
                     receptor::String, weight::Float64)
    post in keys(net.neuron) || begin @error "No target neural population can be found."; exit(1) end
    pre in keys(net.synapse) || begin net.synapse[pre] = Vector{Synapse}() end
    push!(net.synapse[pre], Synapse(pre, post, receptor, weight))
end


function simulate(net::Network, t, dt=0.1)
    # index: each neuron
    net_size = length(net.neuron)
    neu_index = Dict(neu.name=>i for (i, neu) in enumerate(values(net.neuron)))

    # index: each receptor
    receptor_accumulated = accumulate(+, (length(neu.receptor) for neu in values(net.neuron)))
    receptor_index_start = [0; receptor_accumulated[1:end-1]] .+ 1
    receptor_index = Dict{String,Dict{String,Int}}()
    for (i, neu) in enumerate(values(net.neuron))
        receptor_index[neu.name] = Dict{String,Int}()
        for (j, rec) in enumerate(values(neu.receptor))
            receptor_index[neu.name][rec.name] = receptor_index_start[i] + j - 1
        end
    end

    # initialization: v
    v = [neu.rest for neu in values(net.neuron)]  # v: membrane voltage

    # initialization: gating variable (fast synapse)
    s = zeros(receptor_accumulated[end])  # faster when combining with v?

    # initialization: current
    I = zeros(net_size)  # TODO: initialization from event
    I[1] = 0.68  # XXX: temp

    # output:  # TODO: if store_potential=true
    potential = zeros(length(0:dt:t), net_size)
    potential[1, :] .= v

    # simulation
    for (i, _) in enumerate(dt:dt:t)
        for (j, neu) in enumerate(values(net.neuron))

            # update synapse
            receptor_current = 0.0
            for (k, rec) in enumerate(values(neu.receptor))
                receptor_current += -s[receptor_index_start[j] + k-1] * (v[j] - rec.reversal) / 1000

                dsdt = -s[receptor_index_start[j] + k-1] / rec.tau
                s[receptor_index_start[j] + k-1] += dsdt * dt
            end

            # update neuron
            dvdt = -(1 / neu.taum) * (v[j] - neu.rest) + +(receptor_current, I[j]) / neu.Cm  # TODO: use decay factor?
            v[j] += dvdt * dt

            if v[j] >= neu.threshold
                v[j] = neu.reset
                # TODO: store a spike
                # TODO: refractory period

                for syn in net.synapse[neu.name]
                    rec_idx = receptor_index[syn.post][syn.receptor]
                    s[rec_idx] += syn.weight
                end
            end

            potential[i+1, j] = v[j]
        end
    end

    return potential
end
