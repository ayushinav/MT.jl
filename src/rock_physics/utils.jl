@generated function from_nt(::Type{T}, nt::NamedTuple) where T
    fnames = fieldnames(T)
    args = [:(getproperty(nt, $(QuoteNode(f)))) for f in fnames]
    return :(T($(args...)))
end

function from_nt(::Type{Nothing}; nt::NamedTuple)
    (;)
end

function to_nt(s)
    T = typeof(s)
    names = fieldnames(T)
    vals = ntuple(i -> getfield(s, names[i]), length(names))
    NamedTuple{names}(vals)
end

forward(m::Nothing, p; params = (;)) = nothing
