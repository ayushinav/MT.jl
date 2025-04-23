for m in subtypes(AbstractGeophyModel)
    mstring = string(m)
    dstring = mstring*"Distribution"
    eval(Meta.parse("sample_type(d::$dstring) = $mstring"))
    eval(Meta.parse("sample_type(::Type{T}) where {T <:$dstring} = $mstring"))
end

for m in subtypes(AbstractCondModel)
    mstring = string(m)
    dstring = mstring*"Distribution"
    eval(Meta.parse("sample_type(d::$dstring) = $mstring"))
    eval(Meta.parse("sample_type(::Type{T}) where {T <:$dstring} = $mstring"))
end

for m in subtypes(AbstractElasticModel)
    mstring = string(m)
    dstring = mstring*"Distribution"
    eval(Meta.parse("sample_type(d::$dstring) = $mstring"))
    eval(Meta.parse("sample_type(::Type{T}) where {T <:$dstring} = $mstring"))
end

for m in subtypes(AbstractViscousModel)
    mstring = string(m)
    dstring = mstring*"Distribution"
    eval(Meta.parse("sample_type(d::$dstring) = $mstring"))
    eval(Meta.parse("sample_type(::Type{T}) where {T <:$dstring} = $mstring"))
end

for m in subtypes(AbstractAnelasticModel)
    mstring = string(m)
    dstring = mstring*"Distribution"
    eval(Meta.parse("sample_type(d::$dstring) = $mstring"))
    eval(Meta.parse("sample_type(::Type{T}) where {T <:$dstring} = $mstring"))
end

function sample_type(d::model_2phaseDistribution{V, T1, T2, M}) where {V,T1, T2, M}

    if typeof(V) <: Distribution
        v = typeof(rand(d.ϕ))
    else
        v = V
    end
    t1 = sample_type(T1)
    t2 = sample_type(T2)
    m = M
    model_2phase{v, t1, t2, m}
end
# 

to_dist_nt(d::T) where T<: MT.AbstractModelDistribution = (@show "ASA"; to_nt(d))

function to_dist_nt(d::T) where T <: model_2phaseDistribution
    @show "HHHHHHHHHHHHHAAAAAAAAAA"
    m1 = MT.to_nt(d.m1)
    m2 = MT.to_nt(d.m2)
    return (; ϕ = d.ϕ, m1..., m2...)
end

# function to_dist_nt(d::T) where T <: model_multi_rp
#     m1 = MT.to_nt(d.m1)
#     m2 = MT.to_nt(d.m2)
#     return (; ϕ = d.ϕ, m1..., m2...)
# end
