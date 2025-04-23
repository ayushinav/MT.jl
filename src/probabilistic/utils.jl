for m in subtypes(AbstractGeophyModel)
    mstring = string(m)
    dstring = mstring * "Distribution"
    eval(Meta.parse("sample_type(d::$dstring) = $mstring"))
    eval(Meta.parse("sample_type(::Type{T}) where {T <:$dstring} = $mstring"))
end

for m in subtypes(AbstractCondModel)
    mstring = string(m)
    dstring = mstring * "Distribution"
    eval(Meta.parse("sample_type(d::$dstring) = $mstring"))
    eval(Meta.parse("sample_type(::Type{T}) where {T <:$dstring} = $mstring"))
end

for m in subtypes(AbstractElasticModel)
    mstring = string(m)
    dstring = mstring * "Distribution"
    eval(Meta.parse("sample_type(d::$dstring) = $mstring"))
    eval(Meta.parse("sample_type(::Type{T}) where {T <:$dstring} = $mstring"))
end

for m in subtypes(AbstractViscousModel)
    mstring = string(m)
    dstring = mstring * "Distribution"
    eval(Meta.parse("sample_type(d::$dstring) = $mstring"))
    eval(Meta.parse("sample_type(::Type{T}) where {T <:$dstring} = $mstring"))
end

for m in subtypes(AbstractAnelasticModel)
    mstring = string(m)
    dstring = mstring * "Distribution"
    eval(Meta.parse("sample_type(d::$dstring) = $mstring"))
    eval(Meta.parse("sample_type(::Type{T}) where {T <:$dstring} = $mstring"))
end

function sample_type(d::model_2phaseDistribution{V, T1, T2, M}) where {V, T1, T2, M}
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

sample_type(::Type{Nothing}) = Nothing

function sample_type(d::model_multi_rpDistribution{T1, T2, T3, T4}) where {T1, T2, T3, T4}
    model_multi_rp{sample_type(T1), sample_type(T2), sample_type(T3), sample_type(T4)}
end

# NamedTuple manipulation

to_dist_nt(d::T) where {T <: AbstractModelDistribution} = to_nt(d)

function to_dist_nt(d::T) where {T <: model_2phaseDistribution}
    m1 = MT.to_nt(d.m1)
    m2 = MT.to_nt(d.m2)
    return (; ϕ=d.ϕ, m1..., m2...)
end

function to_dist_nt(d::T) where {T <: model_multi_rpDistribution}
    return merge(to_nt(d.cond), to_nt(d.elastic), to_nt(d.visc), to_nt(d.anelastic))
end

to_dist_nt(d::T) where {T <: AbstractResponseDistribution} = to_nt(d)

function to_dist_nt(d::T) where {T <: multi_rp_responseDistribution}
    return merge(to_nt(d.cond), to_nt(d.elastic), to_nt(d.visc), to_nt(d.anelastic))
end
