abstract type AbstractNeuronModel end

abstract type AbstractSolver end
struct Euler <: AbstractSolver end
struct RK4 <: AbstractSolver end
