## response plots

function plot_response(vars, resp::response; errs=zero(resp), plt_type=:lines,
        kwargs...) where {response <: AbstractGeophyResponse}
    k = fieldnames(response)
    f = Figure()
    axs = [Axis(f[i, 1]) for i in eachindex(k)]

    if plt_type === :lines
        for i in eachindex(axs)
            lines!(axs[i], vars, getproperty(resp, k[i]); kwargs...)
            xscale, yscale = get_scales(response, Val{k[i]}())

            axs[i].xscale = xscale
            axs[i].yscale = yscale

            xlabel, ylabel = get_labels(response, Val{k[i]}())
            axs[i].xlabel = xlabel
            axs[i].ylabel = ylabel
        end

    elseif plt_type === :scatter
        for i in eachindex(axs)
            scatter!(axs[i], vars, getproperty(resp, k[i]); kwargs...)
            xscale, yscale = get_scales(response, Val{k[i]}())

            axs[i].xscale = xscale
            axs[i].yscale = yscale

            xlabel, ylabel = get_labels(response, Val{k[i]}())
            axs[i].xlabel = xlabel
            axs[i].ylabel = ylabel
        end
    elseif plt_type === :errors
        for i in eachindex(axs)
            errorbars!(axs[i], vars, getproperty(resp, k[i]),
                getproperty(errs, k[i]) ./ 2; kwargs...)
            xscale, yscale = get_scales(response, Val{k[i]}())

            axs[i].xscale = xscale
            axs[i].yscale = yscale

            xlabel, ylabel = get_labels(response, Val{k[i]}())
            axs[i].xlabel = xlabel
            axs[i].ylabel = ylabel
        end
    end

    return f, axs
end

function plot_response!(axs, vars, resp::response; errs=zero(resp), plt_type=:lines,
        kwargs...) where {response <: AbstractGeophyResponse}
    k = fieldnames(response)
    # axs = [f.content[i] for i in eachindex(k)]

    if plt_type === :lines
        for i in eachindex(axs)
            lines!(axs[i], vars, getproperty(resp, k[i]); kwargs...)
            xscale, yscale = get_scales(response, Val{k[i]}())

            axs[i].xscale = xscale
            axs[i].yscale = yscale
        end

    elseif plt_type === :scatter
        for i in eachindex(axs)
            scatter!(axs[i], vars, getproperty(resp, k[i]); kwargs...)
            xscale, yscale = get_scales(response, Val{k[i]}())

            axs[i].xscale = xscale
            axs[i].yscale = yscale
        end
    elseif plt_type === :errors
        for i in eachindex(axs)
            errorbars!(axs[i], vars, getproperty(resp, k[i]),
                getproperty(errs, k[i]) ./ 2; kwargs...)
            xscale, yscale = get_scales(response, Val{k[i]}())

            axs[i].xscale = xscale
            axs[i].yscale = yscale
        end
    end
end

## model plots

function plot_model(args...; kwargs...)
    fig = Figure()
    ax = Axis(fig[1, 1])

    plot_model!(fig, args...; kwargs...)
    ax.yreversed = true

    fig, ax
end

function plot_model!(ax, model::m_type; half_space_thickness=1.25 * sum(model.h),
        kwargs...) where {m_type}
    # ax = f.content[1]
    m = model.m
    h = model.h

    m_vec = 10 .^ [m[1], m...]
    h_v = cumsum(h)
    h_vec = [1.0f-2, h_v..., half_space_thickness]

    stairs!(ax, m_vec, h_vec; step=:post, kwargs...)
    ax.yreversed = true
    nothing
end

#= 
fig = Figure()
ax = Axis(fig[1,1])
plot_model!(ax, m2)
ax.yscale = log10
ylims!(ax, 1e4, 1e0)
fig
=#
