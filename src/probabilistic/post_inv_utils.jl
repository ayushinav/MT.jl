"""
    get_model_list(chains::chain, mDist::mdist; 
        trans_utils = (m = log_tf, h = lin_tf,)) where {mdist <: AbstractModelDistribution, chain <: Chains}
returns a list of models from the `Chains` variable obtained from [`stochastic_inverse`](@ref)

## Arguments
* `chains` : `Chains` object obtaiend from the `Turing` model
* `mDist` : *a priori* distribution defined before performing stochastic inversion
"""
function get_model_list(chains::chain, mDist::mdist; 
    trans_utils = (m = log_tf, h = lin_tf,)) where {mdist <: AbstractModelDistribution, chain <: Chains}

    model_type = sample(mDist);

    length_vec = zeros(Int64, length(propertynames(mDist)));
    for (i,k) in enumerate(propertynames(mDist))
        length_vec[i] = (typeof(getproperty(mDist, k)) <: Distribution) ? length(rand(getproperty(mDist, k))) : length(getproperty(mDist, k))
    end

    preds= []
    for k in chains.name_map.parameters
        push!(preds, chains[k].data[:])
    end
    pred= hcat(preds...);
    model_list = [];

    for i in 1:size(pred, 1)
        m = [];
        prev_length = 1;
        for (i,k) in enumerate(fieldnames(model_type))
            if typeof(getproperty(mDist, k)) <: Distribution
                next_length = prev_length + length_vec[i]
                push!(m, pred[i, prev_length : next_length -1])
                prev_length = 0 + next_length
            else
                push!(m, getproperty(mDist, k))
            end
        end

        model_sample = model_type(m...);

        for k in fieldnames(model_type)
            if k in keys(trans_utils)
                getfield(model_sample, k) .= broadcast(getproperty(trans_utils[k], :tf), getfield(model_sample, k))
            end
        end

        push!(model_list, model_sample)
    end
    return model_list

end 

"""
    get_ρ_at_z(pred, zs)
returns the values of the model splatted as a vector in `pred` at points defined by `zs`
"""
function get_ρ_at_z(pred, zs)
    n= 1+ length(pred) ÷2;
    h= cumsum(pred[n+1:end]);
    res= zeros(length(zs));
    idx= zs.<= h[1];
    res[idx].= pred[1];
    
    for i in 2:length(h)
        idx = ((h[i-1].< zs .<= h[i]));
        res[idx] .= pred[i];
    end
    idx= zs.> h[end];
    res[idx].= pred[n];
    return res;
end

"""
    pre_image(mDist::mdist, chains; trans_utils = (m = log_tf, h = lin_tf,),
        grid = (m = collect(-1:0.1:5), h = collect(10 .^ (0:0.1:6)))
        ) where {mdist <: AbstractGeophyModelDistribution}
returns the variables required to plot the PDF image of stochastic inversion

## Arguments
- `mDist`: Geophysical model distribution used to obtain the `chain` for `stochastic_inverse`
- `chain`: `Chains` obtained from `stochastic_inverse`

## Keyword Arguments
- `trans_utils`: `NamedTuple` containing `transform_utils` to move to and from computational and model domain
- `grid`: `NamedTuple` containing the points where the stochastic image is evaluated
"""
function pre_image(mDist::mdist, chains::chain; trans_utils = (m = log_tf, h = lin_tf,),
    grid = (m = collect(-1:0.1:5), h = collect(10 .^ (0:0.1:6)))
    ) where {mdist <: AbstractGeophyModelDistribution, chain <: Chains}
    # we know that geophy model will have `m` and `h`

    # the following code is tested only for 1D

    m_length = length(rand(mDist.m));
    h_length = (typeof(mDist.h) <: Distribution) ? length(rand(mDist.h)) : length(mDist.h);

    preds= []
    for k in chains.name_map.parameters
        push!(preds, chains[k].data[:])
    end
    pred= hcat(preds...);

    if size(pred, 2) == m_length
        broadcast!(getproperty(trans_utils[:m], :tf), pred, pred);
    elseif size(pred, 2) == (m_length + h_length)
        broadcast!(getproperty(trans_utils[:m], :tf), view(pred, :, 1:m_length), view(pred, :, 1:m_length));
        broadcast!(getproperty(trans_utils[:h], :tf), view(pred, :, m_length+1:m_length+h_length), view(pred, :, m_length+1:m_length+h_length));
    else
        error("size of parameters from model distribution `mDist` does not match the values from chains.")
    end

    if size(pred, 2) == (m_length)
        return (m = pred, h = [mDist.h..., sum(mDist.h)]), grid, trans_utils, mDist
    else
        m2 = zeros(eltype(h), size(pred, 1), length(grid[:h]))
        for i in 1:size(pred, 1)
            m2[i,:] .= get_ρ_at_z(pred[i,:], grid[:h]);
        end
        return (m = m2, h = grid[:h]), grid, trans_utils, mDist
    end
end


function kde_image(pre_image::NamedTuple, grid::NamedTuple, trans_utils::NamedTuple, mDist::mdist; kwargs...) where {mdist <: AbstractGeophyModelDistribution}

    K(u)= inv(sqrt(2π))*exp(-u^2/2)

    function get_kde(data, xgrid)
        σ= std(data);
        h= 1.06* σ;
        px= zeros(size(xgrid))
        for (i, x) in enumerate(xgrid)
            s= 0;
            for idata in data
                s= s+ K((x- idata)/h);
            end
            px[i]= inv(length(data)* h)* s 
        end
        return px;
    end

    m_length = length(pre_image[:h]);
    kde_preds = zeros(m_length, length(grid[:m]));

    for i in 1:m_length
        kde_preds[i,:].= get_kde(broadcast(getproperty(trans_utils[:m], :itf), pre_image[:m][:, i]), grid[:m]);
        norm_factor= sum(kde_preds[i, :]);
        kde_preds[i,:].= kde_preds[i,:]./norm_factor;
    end
    # return kde_preds

    plt = heatmap(broadcast(getproperty(trans_utils[:m], :tf), grid[:m]),
     cumsum(pre_image[:h]), kde_preds; kwargs...) #, cmap=reverse(cgrad(:grays))) #, clim = (0, 0.1))
    plot!(plt, xlabel = "ρ (Ωm)", ylabel = "depth (m)")
    return plt
end

# WIP

using Statistics
function mean_std_image(pre_image::NamedTuple, grid::NamedTuple, trans_utils::NamedTuple, mDist::mdist; kwargs...) where {mdist <: AbstractGeophyModelDistribution}

    m_length = size(pre_image[:m], 2)
    @show m_length
    μ_m = vec(mean(broadcast(getproperty(trans_utils[:m], :itf), pre_image[:m][:, 1:m_length]), dims = 1))
    σ_m = vec(std(broadcast(getproperty(trans_utils[:m], :itf), pre_image[:m][:, 1:m_length]), dims = 1))

    @show μ_m

    μ₊_m = μ_m .+ σ_m
    μ₋_m = μ_m .- σ_m

    @show length.([pre_image[:h], grid[:h]])

    h_model = copy(pre_image[:h])

    if typeof(mDist.h) <: Distribution # variable discretization
        push!(μ₊_m, μ₊_m[end])
        push!(μ₋_m, μ₋_m[end])
        push!(μ_m, μ_m[end])
    else
        h_model = h_model[1:end-1];
    end
    @show length(h_model), length(μ_m)

    

    μ_model = sample(mDist)(broadcast(getproperty(trans_utils[:m], :tf), μ_m), h_model)
    μ₊_model = sample(mDist)(broadcast(getproperty(trans_utils[:m], :tf), μ₊_m), h_model)
    μ₋_model = sample(mDist)(broadcast(getproperty(trans_utils[:m], :tf), μ₋_m), h_model)
    @show size.([μ_model.m, μ_model.h])

    plt = plot_model(μ_model, color = "blue", label = "μ")
    plot_model!(plt, μ₊_model, color = "green", label = "μ ± 1σ")
    plot_model!(plt, μ₋_model, color = "green", label = false, legend = :outertopright)
    plot!(plt; kwargs...)

    return plt
end
