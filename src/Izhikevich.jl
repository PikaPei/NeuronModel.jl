abstract type NeuronModel end

mutable struct Izhikevich <: NeuronModel
    name::String
    N::Int
    C::Float64
    Vr::Float64
    Vth::Float64
    Vpeak::Float64
    k::Float64
    a::Float64
    b::Float64
    c::Float64
    d::Float64
end

function Izhikevich(name;
                    N=1, C=100.0, Vr=-60.0, Vth=-40.0, Vpeak=35.0,
                    k=0.7, a=0.03, b=-2.0, c=-50.0, d=100.0)
    return Izhikevich(name, N, C, Vr, Vth, Vpeak, k, a, b, c, d)
end
