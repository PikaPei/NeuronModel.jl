abstract type AbstractEvent end

mutable struct EventCurrentInjection <: AbstractEvent
    time::Float64
    type::String
    population::String
    mean::Float64
    std::Float64
end

mutable struct EventExtFreq <: AbstractEvent
    time::Float64
    type::String
    population::String
    receptor::String
    freq::Float64
end

mutable struct EventEndTrial <: AbstractEvent
    time::Float64
    type::String
end
