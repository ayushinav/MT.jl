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

function from_nt(m::Type{T}, ps_nt::NamedTuple) where T<:construct_model_2phase

    # fnames = fieldnames(T)
    ϕ = getproperty(ps_nt, :ϕ)
    m1 = m.types[1].parameters[1]
    m2 = m.types[2].parameters[1]
    mix = m.types[3]

    model1 = MT.from_nt(m1, ps_nt)
    model2 = MT.from_nt(m2, ps_nt)

    return model_2phase(ϕ, model1, model2, mix())

end

forward(m::Nothing, p; params = (;)) = nothing
