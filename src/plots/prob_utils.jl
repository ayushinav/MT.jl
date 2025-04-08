function get_kde_image!(fig,
        chain::C,
        mDist::mdist;
        cb_kwargs=(;),
        hm_kwargs=(;),
        K=gaussian_kernel,
        half_space_thickness=nothing,
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
    kde_img = zeros(length(grid.m), h_length)  # nₘ x nₕ

    for i in 1:h_length
        kde_img[:, i] .= get_kde(pred[:, i], grid.m; K=K)
        norm_factor = sum(kde_img[:, i])
        kde_img[:, i] .= kde_img[:, i] ./ norm_factor
    end

    isnothing(half_space_thickness) && (half_space_thickness = sum(mDist.h) * 1.25)
    zs = [0, cumsum(mDist.h)..., half_space_thickness]

    ax = Axis(fig[1, 1])
    cb = fig[1, 2]

    ms = broadcast(trans_utils.m.tf, grid.m)
    hm = heatmap!(ax, ms, zs, kde_img; hm_kwargs...)
    ax.yreversed = true

    Colorbar(cb, hm, cb_kwargs...)

    if return_kde_mat
        return kde_img
    else
        return nothing
    end
end

function get_kde_image!(fig,
        chain::C,
        mDist::mdist;
        cb_kwargs=(;),
        hm_kwargs=(;),
        K=gaussian_kernel,
        half_space_thickness=nothing,
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

    isnothing(half_space_thickness) && (half_space_thickness = sum(mean(mDist.h)) * 1.25)
    zs = [0, grid.z..., grid.z[end] .+ range(0.0, half_space_thickness; length=10)...]

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
    hm = heatmap!(ax, ms, zs, kde_img; hm_kwargs...)
    ax.yreversed = true

    Colorbar(fig[1, 2], hm, cb_kwargs...)

    if return_kde_mat
        return kde_img
    else
        return nothing
    end
end

function get_kde_image(args...; return_kde_mat=false, kwargs...)
    fig = Figure()
    kde_img = get_kde_image!(fig, args...; kwargs...)

    if return_kde_mat
        return fig, kde_img
    else
        return fig
    end
end

function get_mean_std_image!(ax,
        chain::C,
        mDist::mdist;
        confidence_interval=0.95,
        half_space_thickness=nothing,
        plot_kwargs=nothing,
        trans_utils=(m=lin_tf, h=lin_tf)) where {
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

    isnothing(half_space_thickness) && (half_space_thickness = sum(mDist.h) * 1.25)
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

    m_type = MT.sample(mDist)

    plot_model!(ax, m_type(trans_utils.m.tf.(μ_m), mDist.h); mean_kwargs...)
    plot_model!(ax, m_type(trans_utils.m.tf.(μ₊_m), mDist.h); std_plus_kwargs...)
    plot_model!(ax, m_type(trans_utils.m.tf.(μ₋_m), mDist.h); std_minus_kwargs...)

    nothing
end

function get_mean_std_image!(ax,
        chain::C,
        mDist::mdist;
        confidence_interval=0.95,
        half_space_thickness=nothing,
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

    # broadcast!(getproperty(trans_utils[:m], :tf),
    # view(pred, :, 1:m_length), view(pred, :, 1:m_length))
    broadcast!(trans_utils.h.tf, view(pred, :, (m_length + 1):(m_length + h_length)),
        view(pred, :, (m_length + 1):(m_length + h_length)))

    isnothing(half_space_thickness) && (half_space_thickness = sum(mean(mDist.h)) * 1.25)
    zs = [0, z_points..., z_points[end] .+ range(0.0, half_space_thickness; length=10)...]

    m2 = zeros(eltype(pred), size(pred, 1), length(zs))
    for i in axes(pred, 1)
        m2[i, :] .= MT.get_ρ_at_z(pred[i, :], zs)
    end

    μ_m = mean(m2; dims=1)[:]
    σ_m = mean(m2; dims=1)[:]

    z_score = quantile(Normal(0.0, 1.0), 1 - (1 - confidence_interval) / 2)

    μ₊_m = @. μ_m + z_score * σ_m
    μ₋_m = @. μ_m - z_score * σ_m

    @show μ_m
    @show μ₊_m
    @show μ₋_m

    isnothing(half_space_thickness) && (half_space_thickness = sum(mDist.h) * 1.25)
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

    # plot_model!(ax, m_type(μ_m, mDist.h); mean_kwargs...)
    lines!(ax, trans_utils.m.tf.(μ_m), zs; mean_kwargs...)
    lines!(ax, trans_utils.m.tf.(μ₊_m), zs; mean_kwargs...)
    lines!(ax, trans_utils.m.tf.(μ₋_m), zs; mean_kwargs...)
    # plot_model!(ax, m_type(μ₊_m, mDist.h); std_plus_kwargs...)
    # plot_model!(ax, m_type(μ₋_m, mDist.h); std_minus_kwargs...)

    nothing
end

function get_mean_std_image(args...; kwargs...)
    fig = Figure()
    ax = Axis(fig[1, 1])
    get_mean_std_image!(ax, args...; kwargs...)
    fig
end
