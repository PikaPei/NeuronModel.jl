using Parameters

@with_kw struct Receptor{T<:Number}
    name::String
    tau::T
    reversal::T
end
