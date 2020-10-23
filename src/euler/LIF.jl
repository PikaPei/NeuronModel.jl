using Parameters
include("receptor.jl")

abstract type AbstractNeuronModel end

@with_kw struct LIF{T<:Number} <: AbstractNeuronModel
    name::String
    #index::Int
    N::Int = 1
    Cm::T = 0.5
    taum::T = 15.0
    rest::T = -70.0
    reset::T = -55.0
    threshold::T = -50.0
    receptors::Vector{Receptor} = Vector{Receptor}()
end

LIF(name; kwargs...) = LIF(name=name; kwargs...)

function LIF(name, N, Cm, taum, rest, reset, threshold)
    LIF(name=name, N=N, Cm=Cm, taum=taum,
        rest=rest, reset=reset, threshold=threshold)
end


function simulate(neuron::LIF, u, p, t, dt)
    # p: input current

    potential = zeros(length(0:dt:t))
    potential[1] = neuron.rest

    for (i, _) in enumerate(dt:dt:t)
        dvdt = -(1 / neuron.taum) * (u - neuron.rest) + p / neuron.Cm
        u += dvdt * dt

        if u >= neuron.threshold
            u = neuron.reset
            # fire a spike
        end

        potential[i+1] = u
    end

    return potential
end
