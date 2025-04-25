"""
    gaussian_kernel(u, σ² = 2)

return gaussian kernel centered at 0, at points given by `u` with std deviation of `σ`
"""
gaussian_kernel(u, σ²=2) = inv(sqrt(2π)) * exp(-u^2 / σ²)

"""
    get_kde(data, xgrid; Κ= gaussian_kernel)

returns distribution of `data` using kernel density estimation

## Arguments

  - `data` : 1D vector to evaluate distribution for
  - `xgrid` : 1D vector to evaluate distribution on

## Keyword Arguments

  - `Κ` : kernel used to evaluate density, defaults to `gaussian_kernel`[@ref]
"""
function get_kde(data, xgrid; Κ=gaussian_kernel)
    σ = std(data)
    h = 1.06 * σ
    px = zeros(size(xgrid))
    for (i, x) in enumerate(xgrid)
        s = 0
        for idata in data
            s = s + Κ((x - idata) / h)
        end
        px[i] = inv(length(data) * h) * s
    end
    return px
end

"""
    get_kde_image!(fig,
        chain::C,
        mDist::mdist;
        hm_kwargs=(;),
        cb_kwargs=(;),
        K=gaussian_kernel,
        half_space_depth=nothing,
        kde_transformation_fn = identity,
        return_kde_mat=false,
        trans_utils=(m=lin_tf, h=lin_tf),
        grid=(m=collect(-1:0.1:5), z=cumsum(mDist.h))) where {
        C <: Chains, mdist <: MTModelDistribution{<:Distribution, <:AbstractArray}}

plots on `fig`, a heatmap of probability distributions sampled by a `chain` using kernel density estimation

## Arguments

  - `fig` : Figure on which the heatmap is plotted
  - `chain` : samples in the form `Turing.Chains` from an MCMC sampling
  - `mDist` : *apriori* model distribution used for MCMC sampling

## Keyword Arguments

  - `hm_kwargs` : `NamedTuple` containing keyword arguments for customizing heatmap
  - `cb_kwargs` : `NamedTuple` containing keyword arguments for customizing colorbar
  - `K` : kernel used to perform kernel density estimation
  - `half_space_depth` : extent of half space, i.e., the last layer, informs how far to extend the half space, defaults to `1.25 × last `
  - `kde_transformation_fn` : a function that transforms the image domain, eg., use `log10` to plot log pdf; defaults to `identity` which implies no bounds_transformation
  - `return_kde_mat` : whether to return the matrix containing the values of heatmap along with corresponding x,y axes; defaults to `false`
  - `trans_utils` : `NamedTuple` containing functions to transform the samples; defaults to no `lin_tf` for all parameters
  - `grid` : `NamedTuple` containing grid to evaluate the kernel density on. `m` refers to the points to evaluate kde of model parameters,
    `z` refers to the depth points at which the model samples are inferred, not used if `h` is not sampled.

!!! note


Also check relevant tutorial page!
"""
function get_kde_image!(fig,
        chain::C,
        mDist::mdist;
        cb_kwargs=(;),
        hm_kwargs=(;),
        K=gaussian_kernel,
        half_space_depth=nothing,
        kde_transformation_fn=identity,
        return_kde_mat=false,
        trans_utils=(m=lin_tf, h=lin_tf),
        grid=(m=collect(-1:0.1:5), z=cumsum(mDist.h))) where {
        C <: Chains, mdist <: MTModelDistribution{<:Distribution, <:AbstractArray}}
    preds = []
    for k in chain.name_map.parameters
        push!(preds, chain[k].data[:])
    end
    pred = hcat(preds...)

    h_length = length(mDist.h)
    kde_img = zeros(length(grid.m), h_length + 1)  # nₘ x nₕ

    for i in 1:h_length
        kde_img[:, i] .= get_kde(pred[:, i], grid.m; K=K)
        norm_factor = sum(kde_img[:, i])
        kde_img[:, i] .= kde_img[:, i] ./ norm_factor
    end

    kde_img[:, end] .= copy(kde_img[:, end - 1])

    isnothing(half_space_depth) && (half_space_depth = sum(mDist.h) * 1.25)
    zs = [0, cumsum(mDist.h)..., half_space_depth]

    ax = Axis(fig[1, 1])
    cb = fig[1, 2]

    ms = broadcast(trans_utils.m.tf, grid.m)
    hm = heatmap!(ax, ms, zs, kde_transformation_fn.(kde_img); hm_kwargs...)
    ax.yreversed = true

    Colorbar(cb, hm, cb_kwargs...)

    if return_kde_mat
        return (ms, zs, kde_transformation_fn.(kde_img))
    else
        return nothing
    end
end

"""
    get_kde_image!(fig,
        chain,
        mDist;
        hm_kwargs=(;),
        cb_kwargs=(;),
        K=gaussian_kernel,
        half_space_depth=nothing,
        kde_transformation_fn = identity,
        return_kde_mat=false,
        trans_utils=(m=lin_tf, h=lin_tf),
        grid=(m=collect(-1:0.1:5), z=cumsum(mDist.h))) where {
        C <: Chains, mdist <: MTModelDistribution}

plots on `fig`, a heatmap of probability distributions sampled by a `chain` using kernel density estimation

## Arguments

  - `fig` : Figure on which the heatmap is plotted
  - `chain` : samples in the form `Turing.Chains` from an MCMC sampling
  - `mDist` : *apriori* model distribution used for MCMC sampling

## Keyword Arguments

  - `hm_kwargs` : `NamedTuple` containing keyword arguments for customizing heatmap
  - `cb_kwargs` : `NamedTuple` containing keyword arguments for customizing colorbar
  - `K` : kernel used to perform kernel density estimation
  - `half_space_depth` : extent of half space, i.e., the last layer, informs how far to extend the half space, defaults to `1.25 × last `
  - `kde_transformation_fn` : a function that transforms the image domain, eg., use `log10` to plot log pdf; defaults to `identity` which implies no bounds_transformation
  - `return_kde_mat` : whether to return the matrix containing the values of heatmap along with corresponding x,y axes; defaults to `false`
  - `trans_utils` : `NamedTuple` containing functions to transform the samples; defaults to no `lin_tf` for all parameters
  - `grid` : `NamedTuple` containing grid to evaluate the kernel density on. `m` refers to the points to evaluate kde of model parameters,
    `z` refers to the depth points at which the model samples are inferred, not used if `h` is not sampled.

!!! note


Also check relevant tutorial page!
"""
function get_kde_image!(fig,
        chain::C,
        mDist::mdist;
        cb_kwargs=(;),
        hm_kwargs=(;),
        K=gaussian_kernel,
        half_space_depth=nothing,
        kde_transformation_fn=identity,
        return_kde_mat=false,
        trans_utils=(m=lin_tf, h=lin_tf),
        grid=(m=collect(-1:0.1:5), z=cumsum(mean(mDist.h)))) where {
        C <: Chains, mdist <: MTModelDistribution{<:Distribution, <:Distribution}}
    preds = []
    for k in chain.name_map.parameters
        push!(preds, chain[k].data[:])
    end
    pred = hcat(preds...)

    m_length = length(rand(mDist.m))
    h_length = length(rand(mDist.h))

    broadcast!(trans_utils.m.tf, view(pred, :, 1:m_length), view(pred, :, 1:m_length))
    broadcast!(trans_utils.h.tf, view(pred, :, (m_length + 1):(m_length + h_length)),
        view(pred, :, (m_length + 1):(m_length + h_length)))

    isnothing(half_space_depth) && (half_space_depth = sum(mean(mDist.h)) * 1.25)
    zs = [0, grid.z..., grid.z[end] .+ range(0.0, half_space_depth; length=10)...]

    m2 = zeros(eltype(pred), size(pred, 1), length(zs))
    for i in axes(pred, 1)
        m2[i, :] .= MT.get_ρ_at_z(pred[i, :], zs)
    end

    z_length = length(zs)
    kde_img = zeros(length(grid.m), z_length)  # nₘ x nₕ

    for i in 1:z_length
        kde_img[:, i] .= get_kde(m2[:, i], grid.m; K=K)
        norm_factor = sum(kde_img[:, i])
        kde_img[:, i] .= kde_img[:, i] ./ norm_factor
    end

    ax = Axis(fig[1, 1])
    cb = fig[1, 2]

    ms = broadcast(trans_utils.m.tf, grid.m)
    hm = heatmap!(ax, ms, zs, kde_transformation_fn.(kde_img); hm_kwargs...)
    ax.yreversed = true

    Colorbar(fig[1, 2], hm, cb_kwargs...)

    if return_kde_mat
        return (ms, zs, kde_transformation_fn.(kde_img))
    else
        return nothing
    end
end

"""
    get_kde_image!(fig,
        chain,
        mDist;
        hm_kwargs=(;),
        cb_kwargs=(;),
        K=gaussian_kernel,
        half_space_depth=nothing,
        kde_transformation_fn = identity,
        return_kde_mat=false,
        trans_utils=(m=lin_tf, h=lin_tf),
        grid=(m=collect(-1:0.1:5), z=cumsum(mDist.h)))

returns `fig`, a heatmap of probability distributions sampled by a `chain` using kernel density estimation

## Arguments

  - `fig` : Figure on which the heatmap is plotted
  - `chain` : samples in the form `Turing.Chains` from an MCMC sampling
  - `mDist` : *apriori* model distribution used for MCMC sampling

## Keyword Arguments

  - `hm_kwargs` : `NamedTuple` containing keyword arguments for customizing heatmap
  - `cb_kwargs` : `NamedTuple` containing keyword arguments for customizing colorbar
  - `K` : kernel used to perform kernel density estimation
  - `half_space_depth` : extent of half space, i.e., the last layer, informs how far to extend the half space, defaults to `1.25 × last `
  - `kde_transformation_fn` : a function that transforms the image domain, eg., use `log10` to plot log pdf; defaults to `identity` which implies no bounds_transformation
  - `return_kde_mat` : whether to return the matrix containing the values of heatmap along with corresponding x,y axes; defaults to `false`
  - `trans_utils` : `NamedTuple` containing functions to transform the samples; defaults to no `lin_tf` for all parameters
  - `grid` : `NamedTuple` containing grid to evaluate the kernel density on. `m` refers to the points to evaluate kde of model parameters,
    `z` refers to the depth points at which the model samples are inferred, not used if `h` is not sampled.

!!! note


Also check relevant tutorial page!
"""
function get_kde_image(args...; return_kde_mat=false, kwargs...)
    fig = Figure()
    kde_img = get_kde_image!(fig, args...; kwargs...)

    if return_kde_mat
        return fig, kde_img
    else
        return fig
    end
end

"""
    get_mean_std_image!(ax,
        chain,
        mDist::mdist;
        confidence_interval=0.95,
        half_space_depth=nothing,
        plot_kwargs=nothing,
        trans_utils=(m=lin_tf, h=lin_tf))

plots on `ax`, a bounds plot (using mean and std deviation) of probability distributions sampled by a `chain` using kernel density estimation

## Arguments

  - `fig` : Axis on which the probability bounds are plotted
  - `chain` : samples in the form `Turing.Chains` from an MCMC sampling
  - `mDist` : *apriori* model distribution used for MCMC sampling

## Keyword Arguments

  - `confidence_interval` : a `confidence_interval` of `0.9` implies 90% of values are within the bounds
  - `half_space_depth` : extent of half space, i.e., the last layer, informs how far to extend the half space, defaults to `1.25 × last `
  - `plot_kwargs` : `NamedTuple` containing keyword arguments for plots
  - `return_kde_mat` : whether to return the matrix containing the values of heatmap along with corresponding x,y axes; defaults to `false`
  - `trans_utils` : `NamedTuple` containing functions to transform the samples; defaults to no `lin_tf` for all parameters
  - `z_points` : depth points at which bounds are plotted, not used if `h` is not sampled; defaults to depths corresponding to `mean(h)`

!!! note


Also check relevant tutorial page!
"""
function get_mean_std_image!(ax,
        chain,
        mDist;
        confidence_interval=0.95,
        half_space_depth=nothing,
        plot_kwargs=nothing,
        trans_utils=(m=lin_tf, h=lin_tf),
        z_points=cumsum(mDist.h)) where {
        C <: Chains, mdist <: MTModelDistribution{<:Distribution, <:AbstractArray}}
    preds = []
    for k in chain.name_map.parameters
        push!(preds, chain[k].data[:])
    end
    pred = hcat(preds...)

    μ_m = mean(pred; dims=1)[:]
    σ_m = mean(pred; dims=1)[:]

    z_score = quantile(Normal(0.0, 1.0), 1 - (1 - confidence_interval) / 2)

    μ₊_m = @. μ_m + z_score * σ_m
    μ₋_m = @. μ_m - z_score * σ_m

    isnothing(half_space_depth) && (half_space_depth = sum(mDist.h) * 1.25)
    mean_kwargs = (;)
    std_plus_kwargs = (;)
    std_minus_kwargs = (;)
    if isnothing(plot_kwargs)
        mean_kwargs = (mean_kwargs..., label="mean", color=:blue)

        ci = Int(100confidence_interval)

        std_plus_kwargs = (std_plus_kwargs..., label="upper bound", color=:green)

        std_minus_kwargs = (std_minus_kwargs..., label="lower bound", color=:green)
    else
        mean_kwargs = (; plot_kwargs.mean_kwargs...)
        std_plus_kwargs = (; plot_kwargs.std_plus_kwargs...)
        std_minus_kwargs = (; plot_kwargs.std_minus_kwargs...)
    end

    m_type = MT.sample(mDist)

    plot_model!(ax, m_type(trans_utils.m.tf.(μ_m), mDist.h); mean_kwargs...)
    plot_model!(ax, m_type(trans_utils.m.tf.(μ₊_m), mDist.h); std_plus_kwargs...)
    plot_model!(ax, m_type(trans_utils.m.tf.(μ₋_m), mDist.h); std_minus_kwargs...)

    nothing
end

"""
    get_mean_std_image!(ax,
        chain,
        mDist;
        confidence_interval=0.95,
        half_space_depth=nothing,
        plot_kwargs=nothing,
        trans_utils=(m=lin_tf, h=lin_tf))

plots on `ax`, a bounds plot (using mean and std deviation) of probability distributions sampled by a `chain` using kernel density estimation

## Arguments

  - `fig` : Axis on which the probability bounds are plotted
  - `chain` : samples in the form `Turing.Chains` from an MCMC sampling
  - `mDist` : *apriori* model distribution used for MCMC sampling

## Keyword Arguments

  - `confidence_interval` : a `confidence_interval` of `0.9` implies 90% of values are within the bounds
  - `half_space_depth` : extent of half space, i.e., the last layer, informs how far to extend the half space, defaults to `1.25 × last `
  - `plot_kwargs` : `NamedTuple` containing keyword arguments for plots
  - `return_kde_mat` : whether to return the matrix containing the values of heatmap along with corresponding x,y axes; defaults to `false`
  - `trans_utils` : `NamedTuple` containing functions to transform the samples; defaults to no `lin_tf` for all parameters
  - `z_points` : depth points at which bounds are plotted, not used if `h` is not sampled; defaults to depths corresponding to `mean(h)`

!!! note


Also check relevant tutorial page!
"""
function get_mean_std_image!(ax,
        chain::C,
        mDist::mdist;
        confidence_interval=0.95,
        half_space_depth=nothing,
        plot_kwargs=nothing,
        trans_utils=(m=lin_tf, h=lin_tf),
        z_points=cumsum(mean(mDist.h))) where {
        C <: Chains, mdist <: MTModelDistribution{<:Distribution, <:Distribution}}
    preds = []
    for k in chain.name_map.parameters
        push!(preds, chain[k].data[:])
    end
    pred = hcat(preds...)

    m_length = length(rand(mDist.m))
    h_length = length(rand(mDist.h))

    broadcast!(trans_utils.h.tf, view(pred, :, (m_length + 1):(m_length + h_length)),
        view(pred, :, (m_length + 1):(m_length + h_length)))

    isnothing(half_space_depth) && (half_space_depth = sum(mean(mDist.h)) * 1.25)
    zs = [0, z_points..., z_points[end] .+ range(0.0, half_space_depth; length=10)...]

    m2 = zeros(eltype(pred), size(pred, 1), length(zs))
    for i in axes(pred, 1)
        m2[i, :] .= MT.get_ρ_at_z(pred[i, :], zs)
    end

    μ_m = mean(m2; dims=1)[:]
    σ_m = mean(m2; dims=1)[:]

    z_score = quantile(Normal(0.0, 1.0), 1 - (1 - confidence_interval) / 2)

    μ₊_m = @. μ_m + z_score * σ_m
    μ₋_m = @. μ_m - z_score * σ_m

    isnothing(half_space_depth) && (half_space_depth = sum(mDist.h) * 1.25)
    mean_kwargs = (;)
    std_plus_kwargs = (;)
    std_minus_kwargs = (;)
    if isnothing(plot_kwargs)
        mean_kwargs = (mean_kwargs..., label="mean", color=:blue)

        std_plus_kwargs = (std_plus_kwargs..., label="upper bound", color=:green)

        std_minus_kwargs = (std_minus_kwargs..., label="lower bound", color=:green)
    else
        mean_kwargs = (; plot_kwargs.mean_kwargs...)
        std_plus_kwargs = (; plot_kwargs.std_plus_kwargs...)
        std_minus_kwargs = (; plot_kwargs.std_minus_kwargs...)
    end

    lines!(ax, trans_utils.m.tf.(μ_m), zs; mean_kwargs...)
    lines!(ax, trans_utils.m.tf.(μ₊_m), zs; std_plus_kwargs...)
    lines!(ax, trans_utils.m.tf.(μ₋_m), zs; std_minus_kwargs...)

    nothing
end

"""
    get_mean_std_image!(chain,
        mDist;
        confidence_interval=0.95,
        half_space_depth=nothing,
        plot_kwargs=nothing,
        trans_utils=(m=lin_tf, h=lin_tf))

return `fig`, a figure with a bounds plot (using mean and std deviation) of probability distributions sampled by a `chain` using kernel density estimation

## Arguments

  - `fig` : Axis on which the probability bounds are plotted
  - `chain` : samples in the form `Turing.Chains` from an MCMC sampling
  - `mDist` : *apriori* model distribution used for MCMC sampling

## Keyword Arguments

  - `confidence_interval` : a `confidence_interval` of `0.9` implies 90% of values are within the bounds
  - `half_space_depth` : extent of half space, i.e., the last layer, informs how far to extend the half space, defaults to `1.25 × last `
  - `plot_kwargs` : `NamedTuple` containing keyword arguments for plots
  - `return_kde_mat` : whether to return the matrix containing the values of heatmap along with corresponding x,y axes; defaults to `false`
  - `trans_utils` : `NamedTuple` containing functions to transform the samples; defaults to no `lin_tf` for all parameters
  - `z_points` : depth points at which bounds are plotted, not used if `h` is not sampled; defaults to depths corresponding to `mean(h)`

!!! note


Also check relevant tutorial page!
"""
function get_mean_std_image(args...; kwargs...)
    fig = Figure()
    ax = Axis(fig[1, 1])
    get_mean_std_image!(ax, args...; kwargs...)
    fig
end
