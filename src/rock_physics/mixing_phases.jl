# mutable struct HS{T where T <: AbstractRockphyModel} <: AbstractRockphyModel
#     models::Vector{T}
# end

function get_forward_function(hs_plus::HS_plus)

    property_lists = [];
    [push!(property_lists, [propertynames(m)...]) for m in models_list]

    plist = [];
    for i in eachindex(property_lists)
        push!(plist, property_lists[i]...)
    end

    plist

    prop_name = [];
    model_dicts = []

    for p in plist
        if p ∈ prop_name
            continue
        end
        k = p .∈ propertynames.(models_list)
        push!(prop_name, p)
        push!(model_dicts, models_list[k])

    end

    # remove duplicates
    i =1
    n = length(plist)
    while i <= n
        # push!(idcs, i);
        j = i+1;
        while (j <= n) && (i <= n)
        # for j in i+1:length(a)
        # @show i, j, a[i], a[j]
            if plist[i] == plist[j]
                deleteat!(plist, j)
                n = n- 1
                j=j-1
            end
            j = j+1;
            # @show a
        end
        i = i+1;
    end




    return function HS_plus(;args...)
        # @show length(args), length(plist)

        @assert length(args) == length(plist)

        m_total = 0.;
        for i in eachindex(hs_plus.models)
            # @show fieldnames(hs_plus.models[i])
            # @show args, typeof(args)
            # @show args[fieldnames(hs_plus.models[i])]
            σi = hs_plus.models[i](args[fieldnames(hs_plus.models[i])...]...)

            # @show σi

            m_total = m_total + hs_plus.ϕ[i] * forward(σi);
        end
        return m_total
    end

end