include("receptor.jl")

@with_kw mutable struct LIF{T<:Real} <: AbstractNeuronModel
    name::String
    #index::Int
    N::Int = 1
    Cm::T = 0.5
    taum::T = 15.0
    rest::T = -70.0
    reset::T = -55.0
    threshold::T = -50.0
    refractory_period = 2.0
    refractory_status = false
    recovery_timing = 0.0
    receptor::Vector{Receptor} = Vector{Receptor}()
end

LIF(name; kwargs...) = LIF(name=name; kwargs...)

function LIF(name, N, Cm, taum, rest, reset, threshold;
             refractory_period=2.0, refractory_status=false,
             recovery_timing=0.0)
    LIF(name=name, N=N, Cm=Cm, taum=taum,
        rest=rest, reset=reset, threshold=threshold,
        refractory_period=refractory_period, refractory_status=refractory_status,
        recovery_timing=recovery_timing)
end


function lif_eq(neu::LIF, v, current)
    # TODO: use decay factor?
    return -(1 / neu.taum) * (v - neu.rest) + current / neu.Cm
end


function lif_solve(solver::Euler, neu::LIF, v, current, dt)
    return v + dt * lif_eq(neu, v, current)
end


function lif_solve!(solver::Euler, neu::LIF, v::Vector, neu_idx, current, dt)
    v[neu_idx] += dt * lif_eq(neu, v[neu_idx], current)
end


function lif_solve(solver::RK4, neu::LIF, v, current, dt)
    k1 = dt * lif_eq(neu, v, current)
    k2 = dt * lif_eq(neu, v + 0.5k1, current)
    k3 = dt * lif_eq(neu, v + 0.5k2, current)
    k4 = dt * lif_eq(neu, v + k3, current)
    return v + 1/6 * (k1 + 2k2 + 2k3 + k4)
end


function lif_solve!(solver::RK4, neu::LIF, v::Vector, neu_idx, current, dt)
    _v = v[neu_idx]
    k1 = dt * lif_eq(neu, _v, current)
    k2 = dt * lif_eq(neu, _v + 0.5k1, current)
    k3 = dt * lif_eq(neu, _v + 0.5k2, current)
    k4 = dt * lif_eq(neu, _v + k3, current)
    v[neu_idx] += (1/6 * (k1 + 2k2 + 2k3 + k4))
end


# """
# single neuron test
# """
# function simulate(neuron::LIF, u, p, t, dt)
#     # p: input current
#
#     potential = zeros(length(0:dt:t))
#     potential[1] = neuron.rest
#
#     for (i, _) in enumerate(dt:dt:t)
#         dvdt = -(1 / neuron.taum) * (u - neuron.rest) + p / neuron.Cm
#         u += dvdt * dt
#
#         if u >= neuron.threshold
#             u = neuron.reset
#             # fire a spike
#         end
#
#         potential[i+1] = u
#     end
#
#     return potential
# end
