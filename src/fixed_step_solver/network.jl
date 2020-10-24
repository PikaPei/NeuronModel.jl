using DataStructures
using Parameters


include("types.jl")
include("receptor.jl")
include("synapse.jl")
include("event.jl")
include("output.jl")
include("LIF.jl")
include("Izhikevich.jl")

default_model = LIF


mutable struct Network
    neuron::OrderedDict{String,AbstractNeuronModel}
    receptor::OrderedDict{String,Receptor}
    synapse::OrderedDict{String,Vector{Synapse}}
    event::Vector{AbstractEvent}
    event_end::EventEnd

    Network() = new(OrderedDict{String,AbstractNeuronModel}(),
                    OrderedDict{String,Receptor}(),
                    OrderedDict{String,Vector{Synapse}}(),
                    Vector{AbstractEvent}(),
                    EventEnd(0.0, "End"))
end

function Base.show(io::IO, net::Network)
    println(io, keys(net.neuron))
end


function add_neuron(net::Network, name, model=default_model; kwargs...)
    neuron = model(name; kwargs...)
    net.neuron[name] = neuron
end


function add_receptor(net::Network,
                      name::String,
                      tau::Real,
                      reversal::Real)
    rec = Receptor(name, tau, reversal)
    net.receptor[name] = rec
end


function set_neuron_param_all(net::Network, N::Int, Cm::Real, taum::Real,
                              rest::Real, reset::Real, threshold::Real)

    for neu in values(net.neuron)
        neu.N = N
        neu.Cm = Cm
        neu.taum = taum
        neu.rest = rest
        neu.reset = reset
        neu.threshold = threshold
    end
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


function add_event(net::Network, time, type::String, args...)
    if type == "Current"
        event = EventCurrentInjection(time, type, args...)
        push!(net.event, event)

    elseif type == "End"
        event = EventEnd(time, type, args...)
        net.event_end = event
    end
end


function simulate(net::Network, ::Type{LIF}; solver=Euler(), dt=0.1, store_potential=false, store_spike=false)
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
    I = zeros(net_size)

    # output
    if store_potential
        potential = zeros(length(0:dt:net.event_end.time), net_size)
    end

    if store_spike
        spike = Vector{Pair{String,Float64}}()
        # sizehint!() ?
    end

    # simulation
    event_index = 1  # TODO: sort event
    for (i, t_now) in enumerate(0:dt:net.event_end.time)
        # update current
        while event_index <= length(net.event) && net.event[event_index].time == t_now
            event = net.event[event_index]
            I[ neu_index[event.population] ] = event.mean
            event_index += 1
        end

        for (j, neu) in enumerate(values(net.neuron))  # TODO: @sync @async -> send synapse will affect?

            # update synapse
            receptor_current = 0.0
            for (k, rec) in enumerate(values(neu.receptor))
                receptor_current += -s[receptor_index_start[j] + k-1] * (v[j] - rec.reversal) / 1000

                dsdt = -s[receptor_index_start[j] + k-1] / rec.tau
                s[receptor_index_start[j] + k-1] += dsdt * dt
            end

            # update neuron
            current = receptor_current + I[j]
            v[j] = lif_solve(solver, neu, v[j], current, dt)

            if v[j] >= neu.threshold
                v[j] = neu.reset

                if store_spike
                    push!(spike, neu.name => t_now)
                end

                # TODO: refractory period

                if haskey(net.synapse, neu.name)
                    for syn in net.synapse[neu.name]
                        rec_idx = receptor_index[syn.post][syn.receptor]
                        s[rec_idx] += syn.weight
                    end
                end
            end

            if store_potential
                potential[i, j] = v[j]
            end
        end
    end

    # return output
    if store_potential && store_spike
        output_spike("SpikeALL.dat", spike, neu_index)
        return (potential, spike)
    elseif store_potential
        return potential
    elseif store_spike
        output_spike("SpikeALL.dat", spike, neu_index)
        return spike
    end
end


function simulate(net::Network, ::Type{Izhikevich}; solver=Euler(), dt=0.1, store_potential=false, store_spike=false)
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

    # initialization: v & u
    izh_v = [neu.Vr for neu in values(net.neuron)]
    izh_u = zeros(net_size)

    # initialization: gating variable (fast synapse)
    s = zeros(receptor_accumulated[end])  # faster when combining with izh_v?

    # initialization: current
    I = zeros(net_size)

    # output
    if store_potential
        potential = zeros(length(0:dt:net.event_end.time), net_size)
    end

    if store_spike
        spike = Vector{Pair{String,Float64}}()
        # sizehint!() ?
    end

    # simulation
    event_index = 1  # TODO: sort event
    for (i, t_now) in enumerate(0:dt:net.event_end.time)
        # update current
        while event_index <= length(net.event) && net.event[event_index].time == t_now
            event = net.event[event_index]
            I[ neu_index[event.population] ] = event.mean
            event_index += 1
        end

        for (j, neu) in enumerate(values(net.neuron))  # TODO: @sync @async -> send synapse will affect?

            # update synapse
            receptor_current = 0.0
            for (k, rec) in enumerate(values(neu.receptor))
                receptor_current += -s[receptor_index_start[j] + k-1] * (izh_v[j] - rec.reversal)  # XXX: without / 1000

                dsdt = -s[receptor_index_start[j] + k-1] / rec.tau
                s[receptor_index_start[j] + k-1] += dsdt * dt
            end

            # update neuron
            v, u = izh_v[j], izh_u[j]
            current = receptor_current + I[j]
            izh_v[j], izh_u[j] = izh_solve(solver, neu, v, u, current, dt)

            if izh_v[j] >= neu.Vpeak
                izh_v[j] = neu.c
                izh_u[j] += neu.d

                if store_spike
                    push!(spike, neu.name => t_now)
                end

                # TODO: refractory period

                if haskey(net.synapse, neu.name)
                    for syn in net.synapse[neu.name]
                        rec_idx = receptor_index[syn.post][syn.receptor]
                        s[rec_idx] += syn.weight
                    end
                end
            end

            if store_potential
                potential[i, j] = izh_v[j]
            end
        end
    end

    # return output
    if store_potential && store_spike
        output_spike("SpikeALL.dat", spike, neu_index)
        return (potential, spike)
    elseif store_potential
        return potential
    elseif store_spike
        output_spike("SpikeALL.dat", spike, neu_index)
        return spike
    end
end


simulate(net::Network; kwargs...) = simulate(net, default_model; kwargs...)
