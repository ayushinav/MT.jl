
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
function pre_image(mDist::mdist, chains::chain; trans_utils=(m=log_tf, h=lin_tf),
                   grid=(m=collect(-1:0.1:5), h=collect(10 .^ (0:0.1:6)))) where {
                                                                                  mdist <:
                                                                                  AbstractGeophyModelDistribution,
                                                                                  chain <:
                                                                                  Chains}
    # we know that geophy model will have `m` and `h`

    # the following code is tested only for 1D

    m_length = length(rand(mDist.m))
    h_length = (typeof(mDist.h) <: Distribution) ? length(rand(mDist.h)) : length(mDist.h)

    preds = []
    for k in chains.name_map.parameters
        push!(preds, chains[k].data[:])
    end
    pred = hcat(preds...)

    if size(pred, 2) == m_length
        broadcast!(getproperty(trans_utils[:m], :tf), pred, pred)
    elseif size(pred, 2) == (m_length + h_length)
        broadcast!(getproperty(trans_utils[:m], :tf), view(pred, :, 1:m_length),
                   view(pred, :, 1:m_length))
        broadcast!(getproperty(trans_utils[:h], :tf),
                   view(pred, :, (m_length + 1):(m_length + h_length)),
                   view(pred, :, (m_length + 1):(m_length + h_length)))
    else
        error("size of parameters from model distribution `mDist` does not match the values from chains.")
    end

    if size(pred, 2) == (m_length)
        return (m=pred, h=[mDist.h..., sum(mDist.h)]), grid, trans_utils, mDist
    else
        m2 = zeros(eltype(grid[:h]), size(pred, 1), length(grid[:h]))
        for i in 1:size(pred, 1)
            m2[i, :] .= get_ρ_at_z(pred[i, :], grid[:h])
        end
        return (m=m2, h=grid[:h]), grid, trans_utils, mDist
    end
end

function get_kde_image(pre_image::NamedTuple, grid::NamedTuple, trans_utils::NamedTuple,
                       mDist::mdist, return_vals=false;
                       kwargs...) where {mdist <: AbstractGeophyModelDistribution}
    K(u) = inv(sqrt(2π)) * exp(-u^2 / 2)

    function get_kde(data, xgrid)
        σ = std(data)
        h = 1.06 * σ
        px = zeros(size(xgrid))
        for (i, x) in enumerate(xgrid)
            s = 0
            for idata in data
                s = s + K((x - idata) / h)
            end
            px[i] = inv(length(data) * h) * s
        end
        return px
    end

    m_length = length(pre_image[:h])
    kde_img = zeros(m_length, length(grid[:m]))

    for i in 1:m_length
        kde_img[i, :] .= get_kde(broadcast(getproperty(trans_utils[:m], :itf),
                                           pre_image[:m][:, i]), grid[:m])
        norm_factor = sum(kde_img[i, :])
        kde_img[i, :] .= kde_img[i, :] ./ norm_factor
    end
    # return kde_preds

    plt = heatmap(broadcast(getproperty(trans_utils[:m], :tf), grid[:m]),
                  cumsum(pre_image[:h]), kde_img; kwargs...) #, cmap=reverse(cgrad(:grays))) #, clim = (0, 0.1))
    plot!(plt; xlabel="ρ (Ωm)", ylabel="depth (m)")
    if return_vals
        return (; kde_img, plt)
    else
        return plt
    end
end

function get_mean_std_image(pre_image::NamedTuple, grid::NamedTuple,
                            trans_utils::NamedTuple, mDist::mdist, return_vals=false;
                            kwargs...) where {mdist <: AbstractGeophyModelDistribution}
    m_length = size(pre_image[:m], 2)
    μ_m = vec(mean(broadcast(getproperty(trans_utils[:m], :itf),
                             pre_image[:m][:, 1:m_length]); dims=1))
    σ_m = vec(std(broadcast(getproperty(trans_utils[:m], :itf),
                            pre_image[:m][:, 1:m_length]); dims=1))

    μ₊_m = μ_m .+ σ_m
    μ₋_m = μ_m .- σ_m

    h_model = copy(pre_image[:h])

    if typeof(mDist.h) <: Distribution # variable discretization
        push!(μ₊_m, μ₊_m[end])
        push!(μ₋_m, μ₋_m[end])
        push!(μ_m, μ_m[end])
    else
        h_model = h_model[1:(end - 1)]
    end

    μ_model = sample(mDist)(broadcast(getproperty(trans_utils[:m], :tf), μ_m), broadcast(getproperty(trans_utils[:h], :tf), h_model))
    μ₊_model = sample(mDist)(broadcast(getproperty(trans_utils[:m], :tf), μ₊_m), broadcast(getproperty(trans_utils[:h], :tf), h_model))
    μ₋_model = sample(mDist)(broadcast(getproperty(trans_utils[:m], :tf), μ₋_m), broadcast(getproperty(trans_utils[:h], :tf), h_model))

    # @show  broadcast(getproperty(trans_utils[:h], :tf), h_model)

    plt = plot_model(μ_model; color="blue", label="μ")
    plot_model!(plt, μ₊_model; color="green", label="μ ± 1σ")
    plot_model!(plt, μ₋_model; color="green", label=false, legend=:outertopright)
    plot!(plt; kwargs...)

    if return_vals
        return (; μ_model, μ₊_model, μ₋_model, plt)
    else
        return plt
    end
end
