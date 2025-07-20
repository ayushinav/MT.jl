for m in subtypes(AbstractGeophyModel)
    mstring = string(m)
    dstring = mstring * "Distribution"
    eval(Meta.parse("sample_type(d::$dstring) = $mstring"))
    eval(Meta.parse("sample_type(::Type{T}) where {T <:$dstring} = $mstring"))
end

for m in subtypes(AbstractCondModel)
    mstring = string(m)
    dstring = mstring * "Distribution"
    eval(Meta.parse("sample_type(d::MT.$dstring) = $mstring"))
    eval(Meta.parse("sample_type(::Type{T}) where {T <:MT.$dstring} = $mstring"))
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

function sample_type(d::two_phase_modelDistribution{V, T1, T2, M}) where {V, T1, T2, M}
    if typeof(V) <: Distribution
        v = typeof(rand(d.ϕ))
    else
        v = V
    end
    t1 = sample_type(T1)
    t2 = sample_type(T2)
    m = M
    two_phase_model{v, t1, t2, m}
end

function sample_type(::Type{two_phase_modelDistribution{V, T1, T2, M}}) where {V, T1, T2, M}
    v = Vector{Float64}
    t1 = sample_type(T1)
    t2 = sample_type(T2)
    m = M
    two_phase_model{v, t1, t2, m}
end

function sample_type(::Type{two_phase_modelDistributionType{T1, T2, M}}) where {T1, T2, M}
    v = Vector{Float64}
    t1 = sample_type(T1)
    t2 = sample_type(T2)
    m = M
    two_phase_model{v, t1, t2, m}
end

sample_type(::Type{Nothing}) = Nothing

function sample_type(::multi_rp_modelDistribution{T1, T2, T3, T4}) where {T1, T2, T3, T4}
    multi_rp_model{sample_type(T1), sample_type(T2), sample_type(T3), sample_type(T4)}
end

function sample_type(::Type{multi_rp_modelDistributionType{
        T1, T2, T3, T4}}) where {T1, T2, T3, T4}
    multi_rp_model{sample_type(T1), sample_type(T2), sample_type(T3), sample_type(T4)}
end

function sample_type(d::tune_rp_modelDistribution{K, M}) where {K, M}
    m_ = sample_type(d.model)
    tune_rp_modelType{Vector{Function}, m_}
end

# NamedTuple manipulation

to_dist_nt(d::T) where {T <: AbstractModelDistribution} = to_nt(d)
to_dist_nt(::Nothing) = (;)
to_nt(::Nothing) = (;)

function to_dist_nt(d::T) where {T <: two_phase_modelDistribution}
    m1 = to_nt(d.m1)
    m2 = to_nt(d.m2)
    return (; ϕ=d.ϕ, m1..., m2...)
end

function to_dist_nt(d::T) where {T <: multi_rp_modelDistribution}
    return merge(to_dist_nt(d.cond), to_dist_nt(d.elastic),
        to_dist_nt(d.visc), to_dist_nt(d.anelastic))
end

to_dist_nt(d::T) where {T <: AbstractResponseDistribution} = to_nt(d)

function to_dist_nt(d::T) where {T <: multi_rp_responseDistribution}
    return merge(to_nt(d.cond), to_nt(d.elastic), to_nt(d.visc), to_nt(d.anelastic))
end

function to_dist_nt(d::tune_rp_modelDistribution)
    return (; d.ps_nt..., fn_list=d.fn_list)
end
