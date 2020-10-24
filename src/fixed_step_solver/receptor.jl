@with_kw struct Receptor{T<:Number}
    name::String
    tau::T
    reversal::T
end


function receptor_eq(r::Receptor, s)
    return -s / r.tau
end


function receptor_solve(solver::Euler, r::Receptor, s, dt)
    return s + dt * receptor_eq(r, s)
end


function receptor_solve(solver::RK4, r::Receptor, s, dt)
    k1 = dt * receptor_eq(r, s)
    k2 = dt * receptor_eq(r, s + 0.5k1)
    k3 = dt * receptor_eq(r, s + 0.5k2)
    k4 = dt * receptor_eq(r, s + k3)
    return s + 1/6 * (k1 + 2k2 + 2k3 + k4)
end


# TODO: use default_solver in the future
receptor_solve(r::Receptor, s, dt) = receptor_solve(Euler(), r, s, dt)
