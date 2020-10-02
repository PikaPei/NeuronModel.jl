using DifferentialEquations
using ModelingToolkit
using DataStructures
using Random
# Random.seed!(1234)

include("network.jl")


@parameters t
@derivatives D'~t


spikes = Vector{Pair{String,Float64}}()
# sizehint!() ?


mutable struct NetworkIndex
    vars::OrderedDict{String,Operation}
    ode_idx::OrderedDict{String,Int}
    param_idx::OrderedDict{String,Int}

    NetworkIndex() = new(OrderedDict{String,Operation}(),
                         OrderedDict{String,Int}(),
                         OrderedDict{String,Int}())
end

network_idx = NetworkIndex()


function gen_variables(net::Network)
    vars = OrderedDict{String,Operation}()
    param_idx = OrderedDict{String,Int}()

    counter = 1
    for neu in values(net.neu)
        v = neu.name * "_v"
        u = neu.name * "_u"
        vars[v] = Variable(Symbol(v))(t)
        vars[u] = Variable(Symbol(u))(t)

        current = neu.name * "_I"
        vars[current] = Variable(Symbol(current))()
        param_idx[current] = counter
        counter += 1

        for rec in neu.receptors
            neu_rec = neu.name * '_' * rec.name
            vars[neu_rec] = Variable(Symbol(neu_rec))(t)
        end
    end

    network_idx.vars = vars
    network_idx.param_idx = param_idx
end


function Izh_v(neu::NeuronModel)
    v = network_idx.vars[neu.name*"_v"]
    u = network_idx.vars[neu.name*"_u"]
    neu.k * (v - neu.Vr) * (v - neu.Vth) - u
end


function Izh_u(neu::NeuronModel)
    v = network_idx.vars[neu.name*"_v"]
    u = network_idx.vars[neu.name*"_u"]
    neu.a * (neu.b * (v - neu.Vr) - u)
end


function Izh_receptor(neu::NeuronModel, rec::Receptor)
    neu_rec = neu.name * '_' * rec.name
    -network_idx.vars[neu_rec] * (network_idx.vars[neu.name*"_v"] - rec.reversal)  # * gReceptor
    # NOTE: unit
end


function gen_ode(net::Network)
    ode_idx = OrderedDict{String,Int}()

    eqs = Equation[]
    counter = 1
    for neu in values(net.neu)

        # receptor
        receptor_current = Operation[]
        for rec in neu.receptors
            neu_rec = neu.name * '_' * rec.name
            eq = D(network_idx.vars[neu_rec]) ~ -network_idx.vars[neu_rec] / rec.tau

            push!(eqs, eq)
            ode_idx[neu_rec] = counter
            counter += 1

            push!(receptor_current, Izh_receptor(neu, rec))
        end

        # neuron
        v_name = neu.name*"_v"
        u_name = neu.name*"_u"
        v = network_idx.vars[v_name]
        u = network_idx.vars[u_name]
        current = network_idx.vars[neu.name*"_I"]

        eq_v = D(v) ~ +(Izh_v(neu), receptor_current..., current) / neu.C
        eq_u = D(u) ~ Izh_u(neu)

        push!(eqs, eq_v)
        ode_idx[v_name] = counter
        counter += 1

        push!(eqs, eq_u)
        ode_idx[u_name] = counter
        counter += 1

    end

    network_idx.ode_idx = ode_idx
    return eqs
end


function callback_fire(neu::NeuronModel)
    condition(u, t, integrator) = (u[ network_idx.ode_idx[neu.name*"_v"] ] - neu.Vpeak)

    function affect!(integrator)
        push!(spikes, neu.name => integrator.t)  # take off round()  # TODO: need a better data structure

        integrator.u[ network_idx.ode_idx[neu.name*"_v"] ] = neu.c
        integrator.u[ network_idx.ode_idx[neu.name*"_u"] ] += neu.d

        for tar in values(neu.targets)
            integrator.u[ network_idx.ode_idx[tar.name*'_'*tar.target_receptor] ] += (tar.mean_efficacy * tar.weight)  # No upper bond
        end
    end

    ContinuousCallback(condition, affect!, save_positions=(false,false))
end


function callback_event(event::EventCurrentInjection)
    condition(u, t, integrator) = (t == event.time)

    affect!(integrator) =
        integrator.p[ network_idx.param_idx[event.population*"_I"] ] = (event.Gauss_mean + event.Gauss_std * randn())

    DiscreteCallback(condition, affect!, save_positions=(false,false))
end


function initialize(var::String, net::Network)
    endswith(var, "_v") ? net.neu[var[1:end-2]].Vr : 0.0
end


function gen_tstops(net::Network)
    unique(event.time for event in net.event) |> sort
end


function gen_problem(net::Network)
    gen_variables(net)
    eqs = gen_ode(net)
    de = ODESystem(eqs)
    f = ODEFunction(de,
                    [network_idx.vars[ode] for ode in keys(network_idx.ode_idx)],
                    [network_idx.vars[p] for p in keys(network_idx.param_idx)])

    net.event_end.time != 0.0 || begin @error "You haven't set the end time yet!"; exit(1) end
    # frequency_to_spike(net)

    # Callback Functions
    event_callbacks = [callback_event(event) for event in net.event]
    fire_callbacks = [callback_fire(neu) for neu in values(net.neu)]
    cb = CallbackSet(event_callbacks..., fire_callbacks...)

    u0 = [initialize(ode, net) for ode in keys(network_idx.ode_idx)]
    tspan = (0.0, net.event_end.time)
    p = repeat([0.0], length(net.neu))

    prob = ODEProblem(f, u0, tspan, p, callback=cb)
    save_idxs = [network_idx.ode_idx[neu.name*"_v"] for neu in values(net.neu)]
    (prob, gen_tstops(net), save_idxs)

    # # MethodError
    # sol = solve(prob, tstops=tstops)
    # sol
end


function gen_neuron_index(net)
    idx = Dict{String,Int}()
    for (index, neu) in enumerate(values(net.neu))
        idx[neu.name] = index
    end
    return idx
end


function output_mempot(filename::String, net, sol; dt=0.0, header=true)
    # TODO: output whole table directly
    f = open(filename, "w")
    header == true && println(f, "time", " ", join((neu.name for neu in values(net.neu)), " "))
    if dt == 0.0
        sol_round = [ [round(mp, digits=3) for mp in u] for u in sol.u ]
        for (index, time) in enumerate(sol.t)
            println(f, time, " ", join(sol_round[index], " "))
        end
    else
        for time in 0.0:dt:net.event_end.time
            sol_round = [round(mp, digits=3) for mp in sol(time)]
            println(f, time, " ", join(sol_round, " "))
        end
    end
    close(f)
end


function output_spike(filename::String, net; name=false, digits=3)
    idx = gen_neuron_index(net)
    f = open(filename, "w")
    spikes_round = (spike.first => round(spike.second, digits=digits) for spike in spikes)
    for spike in spikes_round
        if name == true
            println(f, spike.second, " ", spike.first, " ", idx[spike.first])
        else
            println(f, spike.second, " ", idx[spike.first])
        end
    end
    close(f)
end
