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


function receptor_current(::Type{LIF}, r::Receptor, s::Vector, rec_idx, v::Vector, neu_idx)
    -s[rec_idx] * (v[neu_idx] - r.reversal) / 1000
end


function receptor_current(::Type{Izhikevich}, r::Receptor, s::Vector, rec_idx, v::Vector, neu_idx)
    -s[rec_idx] * (v[neu_idx] - r.reversal)  # XXX: without / 1000
end


function receptor_current(::Type{LIF}, r::Receptor, s, v)
    -s * (v - r.reversal) / 1000
end


function receptor_current(::Type{Izhikevich}, r::Receptor, s, v)
    -s * (v - r.reversal)  # XXX: without / 1000
end


function send_spike(net, pre, s, receptor_index)
    for syn in net.synapse[pre]
        rec_idx = receptor_index[syn.post][syn.receptor]
        s[rec_idx] += syn.weight
    end
end


function simulate(net::Network, model::Type{LIF}; solver=Euler(), dt=0.1, store_potential=false, store_spike=false)
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
    s = zeros(receptor_accumulated[end])  # s: gating variable (fast synapse)
    s_temp = zeros(receptor_accumulated[end])
    I = zeros(net_size)  # I: current

    # output
    if store_potential
        potential = zeros(length(0:dt:net.event_end.time), net_size)
    end

    if store_spike
        spike = Vector{Pair{String,Float64}}()
        # sizehint!() ?
    end

    # simulation
    t_steps = 0.0:dt:net.event_end.time
    sort!(net.event, by = x -> x.time)
    event_index = 1
    for (i, t_now) in enumerate(t_steps)
        # update current
        while event_index <= length(net.event) && net.event[event_index].time == t_now
            event = net.event[event_index]
            I[ neu_index[event.population] ] = event.mean
            event_index += 1
        end

        for (j, neu) in enumerate(values(net.neuron))  # TODO: @sync @async -> send synapse will affect?

            # update synapse
            rec_curr = 0.0
            for (k, rec) in enumerate(values(neu.receptor))
                rec_idx = receptor_index_start[j] + k-1

                # rec_curr += receptor_current(model, rec, s[rec_idx], v[j])
                # s[rec_idx] = receptor_solve(rec, s[rec_idx], dt)

                rec_curr += receptor_current(model, rec, s, rec_idx, v, j)
                receptor_solve!(rec, s, rec_idx, dt)
            end

            if (neu.refractory_status == true) && (t_now >= neu.recovery_timing)
                neu.refractory_status = false
            elseif neu.refractory_status == true
                continue
            end

            # update neuron
            current = rec_curr + I[j]
            # v[j] = lif_solve(solver, neu, v[j], current, dt)
            lif_solve!(solver, neu, v, j, current, dt)

            if v[j] >= neu.threshold
                v[j] = neu.reset

                if store_spike
                    push!(spike, neu.name => t_now)
                end

                neu.refractory_status = true
                neu.recovery_timing = t_now + neu.refractory_period

                if haskey(net.synapse, neu.name)
                    send_spike(net, neu.name, s_temp, receptor_index)
                end
            end
        end

        s += s_temp
        s_temp .= 0

        if store_potential
            potential[i, :] .= v
        end
    end

    # return output  # TODO: output_potential
    if store_potential && store_spike
        potential = [t_steps potential]
        output_potential("MemPotALL.dat", potential)
        output_spike("SpikeALL.dat", spike, neu_index)
        return (potential, spike)
    elseif store_potential
        potential = [t_steps potential]
        output_potential("MemPotALL.dat", potential)
        return potential
    elseif store_spike
        output_spike("SpikeALL.dat", spike, neu_index)
        return spike
    end
end


function simulate(net::Network, model::Type{Izhikevich}; solver=Euler(), dt=0.1, store_potential=false, store_spike=false)
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
    s_temp = zeros(receptor_accumulated[end])

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
    t_steps = 0.0:dt:net.event_end.time
    sort!(net.event, by = x -> x.time)
    event_index = 1
    for (i, t_now) in enumerate(t_steps)
        # update current
        while event_index <= length(net.event) && net.event[event_index].time == t_now
            event = net.event[event_index]
            I[ neu_index[event.population] ] = event.mean
            event_index += 1
        end

        for (j, neu) in enumerate(values(net.neuron))  # TODO: @sync @async -> send synapse will affect?

            # update synapse
            rec_curr = 0.0
            for (k, rec) in enumerate(values(neu.receptor))
                rec_idx = receptor_index_start[j] + k-1

                # rec_curr += receptor_current(model, rec, s[rec_idx], izh_v[j])
                # s[rec_idx] = receptor_solve(rec, s[rec_idx], dt)

                rec_curr += receptor_current(model, rec, s, rec_idx, izh_v, j)
                receptor_solve!(rec, s, rec_idx, dt)
            end

            # update neuron
            # v, u = izh_v[j], izh_u[j]
            current = rec_curr + I[j]
            # izh_v[j], izh_u[j] = izh_solve(solver, neu, v, u, current, dt)
            izh_solve!(solver, neu, izh_v, izh_u, j, current, dt)

            if izh_v[j] >= neu.Vpeak
                izh_v[j] = neu.c
                izh_u[j] += neu.d

                if store_spike
                    push!(spike, neu.name => t_now)
                end

                if haskey(net.synapse, neu.name)
                    send_spike(net, neu.name, s_temp, receptor_index)
                end
            end
        end

        s += s_temp
        s_temp .= 0

        if store_potential
            potential[i, :] .= izh_v
        end
    end

    # return output  # TODO: output_potential
    if store_potential && store_spike
        potential = [t_steps potential]
        output_potential("MemPotALL.dat", potential)
        output_spike("SpikeALL.dat", spike, neu_index)
        return (potential, spike)
    elseif store_potential
        potential = [t_steps potential]
        output_potential("MemPotALL.dat", potential)
        return potential
    elseif store_spike
        output_spike("SpikeALL.dat", spike, neu_index)
        return spike
    end
end


simulate(net::Network; kwargs...) = simulate(net, default_model; kwargs...)
