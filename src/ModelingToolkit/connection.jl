"""
Receptor
"""
mutable struct Receptor
    name::String
    tau::Float64
    reversal::Float64
    external_freq::Float64
    mean_external_efficacy::Float64
    mean_external_connection::Float64
end


"""
Target
"""
mutable struct Target
    name::String
    target_receptor::String
    mean_efficacy::Float64
    weight::Float64
end
