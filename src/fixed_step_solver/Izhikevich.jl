@with_kw mutable struct Izhikevich{T<:Real} <: AbstractNeuronModel
    name::String
    # index::Int
    N::Int = 1
    C::T = 100.0
    Vr::T = -60.0
    Vth::T = -40.0
    Vpeak::T = 35.0
    k::T = 0.7
    a::T = 0.03
    b::T = -2.0
    c::T = -50.0
    d::T = 100.0
    receptor::Vector{Receptor} = Vector{Receptor}()
end

Izhikevich(name; kwargs...) = Izhikevich(name=name; kwargs...)


function izh_eq_v(neu::Izhikevich, v, u, current)
    return (neu.k * (v - neu.Vr) * (v - neu.Vth) - u + current) / neu.C
end


function izh_eq_u(neu::Izhikevich, v, u)
    return neu.a * (neu.b * (v - neu.Vr) - u)
end


function izh_solve(solver::Euler, neu::Izhikevich, v, u, current, dt)
    v_next = v + dt * izh_eq_v(neu, v, u, current)
    u_next = u + dt * izh_eq_u(neu, v, u)
    return (v_next, u_next)
end


function izh_solve!(solver::Euler, neu::Izhikevich, izh_v::Vector, izh_u::Vector, neu_idx, current, dt)
    v, u = izh_v[neu_idx], izh_u[neu_idx]
    izh_v[neu_idx] += dt * izh_eq_v(neu, v, u, current)
    izh_u[neu_idx] += dt * izh_eq_u(neu, v, u)
end


function izh_solve(solver::RK4, neu::Izhikevich, v, u, current, dt)
    k1_v = dt * izh_eq_v(neu, v, u, current)
    k1_u = dt * izh_eq_u(neu, v, u)
    k2_v = dt * izh_eq_v(neu, v + 0.5k1_v, u + 0.5k1_u, current)
    k2_u = dt * izh_eq_u(neu, v + 0.5k1_v, u + 0.5k1_u)
    k3_v = dt * izh_eq_v(neu, v + 0.5k2_v, u + 0.5k2_u, current)
    k3_u = dt * izh_eq_u(neu, v + 0.5k2_v, u + 0.5k2_u)
    k4_v = dt * izh_eq_v(neu, v + k3_v, u + k3_u, current)
    k4_u = dt * izh_eq_u(neu, v + k3_v, u + k3_u)
    v_next = v + 1/6 * (k1_v + 2k2_v + 2k3_v + k4_v)
    u_next = u + 1/6 * (k1_u + 2k2_u + 2k3_u + k4_u)
    return (v_next, u_next)
end
