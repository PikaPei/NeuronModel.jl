using Parameters
include("receptor.jl")

abstract type AbstractNeuronModel end

@with_kw struct Izhikevich{T<:Number} <: AbstractNeuronModel
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
    receptors::Vector{Receptor} = Vector{Receptor}()
end

Izhikevich(name; kwargs...) = Izhikevich(name=name; kwargs...)
