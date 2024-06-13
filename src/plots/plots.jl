# most geophysical surveys generally have the same units across, what would be the point of passing units then?

"""
    prepare_plot(d::resp1, ω; plot_type = :scatter, field_names::Vector{Symbol}= [k for k ∈ fieldnames(typeof(d))], 
    kwargs...) where {resp1 <: AbstractResponse}
Prepares and returns a vector of plots for a response `d` (y-axis) parametrized on `ω` (x-axis)

## Arguments
* `d` : The `response` to be plotted 
* `ω` : `d` is usually plotted as a function of some `ω`

## Keyword arguments
* `plot_type` : `:plot` for line plot, `:scatter` for scatter plot. 
* `field_names` : If the response contains more than one parameters, say phase and amplitude, `fieldn_names` determines which of them should be plotted. Defaults to plotting all of them. 
* `kwargs` : Controls properties of the plot, splats the argument to `Plots.jl`
"""
function prepare_plot(d::resp1, ω; plot_type = :scatter, field_names::Vector{Symbol}= [k for k ∈ fieldnames(typeof(d))], 
    kwargs...) where {resp1 <: AbstractResponse}
 
    units= get_units(d);
    yscale= get_scale(d);
    plts= [];
 
    for i in eachindex(field_names)
        pi= getfield(Plots, plot_type)(2π./ω, getfield(d, field_names[i]), xscale=:log10, 
            # yerr = getfield(d_err, field_names[i]),
            yscale= yscale[i], 
            ylabel= "$(units[i])",
            xlabel= "T (s)"; kwargs...
        )
        push!(plts, pi);
    end
 
    return plts
 end

"""
 prepare_plot!(plts::Vector{Any}, d::resp1, ω; plot_type = :scatter, field_names::Vector{Symbol}= [k for k ∈ fieldnames(typeof(d))], 
 kwargs...) where {resp1 <: AbstractResponse}
Modifies `plts`, the vector of plots, with same functionality as `prepare_plot(d, ω)`.
"""
function prepare_plot!(plts::Vector{Any}, d::resp1, ω; plot_type = :scatter, field_names::Vector{Symbol}= [k for k ∈ fieldnames(typeof(d))],
    kwargs...) where {resp1 <: AbstractResponse}

   units= get_units(d);
   yscale= get_scale(d);
   plot_type = Symbol("$(plot_type)!");

   for i in eachindex(field_names)
       Plots.plot!(plts[i], 2π./ω, getfield(d, field_names[i]), xscale=:log10, 
           # yerr = getfield(d_err, field_names[i]),
           yscale= yscale[i], 
           ylabel= "$(units[i])",
           xlabel= "T (s)"; kwargs...
       )
   end
   nothing;
end

"""
    prepare_plot(d::resp1, ω, d_err::resp2; plot_type = :scatter, field_names::Vector{Symbol}= [k for k ∈ fieldnames(typeof(d))], 
    kwargs...) where {resp1 <: AbstractResponse, resp2 <: AbstractResponse}
Prepares and returns a vector of plots for a response `d` (y-axis) parametrized on `ω` (x-axis)

## Arguments
* `d` : The `response` to be plotted 
* `ω` : `d` is usually plotted as a function of some `ω` 
`d_err` : Of the same type as `d`, contains the error bars in response 

## Keyword arguments
* `plot_type` : `:plot` for line plot, `:scatter` for scatter plot. 
* `field_names` : If the response contains more than one parameters, say phase and amplitude, `fieldn_names` determines which of them should be plotted. Defaults to plotting all of them. 
* `kwargs` : Controls properties of the plot, splats the argument to `Plots.jl` 
"""
function prepare_plot(d::resp1, ω, d_err::resp2; plot_type = :scatter, field_names::Vector{Symbol}= [k for k ∈ fieldnames(typeof(d))], 
    kwargs...) where {resp1 <: AbstractResponse, resp2 <: AbstractResponse}

    units= get_units(d);
    yscale= get_scale(d);
    plts= [];

    for i in eachindex(field_names)
        pi= getfield(Plots, plot_type)(2π./ω, getfield(d, field_names[i]), xscale=:log10, 
            yerr = getfield(d_err, field_names[i]),
            yscale= yscale[i], 
            ylabel= "$(units[i])",
            xlabel= "T (s)"; kwargs...
        )
        push!(plts, pi);
    end

    return plts
end

"""
prepare_plot!(plts::Vector{Any}, d::resp1, ω, d_err::resp2; plot_type = :scatter, field_names::Vector{Symbol}= [k for k ∈ fieldnames(typeof(d))],
kwargs...) where {resp1 <: AbstractResponse, resp2 <: AbstractResponse}
Modifies `plts`, the vector of plots, with same functionality as `prepare_plot(d, ω, d_err)`.
"""
function prepare_plot!(plts::Vector{Any}, d::resp1, ω, d_err::resp2; plot_type = :scatter, field_names::Vector{Symbol}= [k for k ∈ fieldnames(typeof(d))],
     kwargs...) where {resp1 <: AbstractResponse, resp2 <: AbstractResponse}

    units= get_units(d);
    yscale= get_scale(d);
    plot_type = Symbol("$(plot_type)!");

    for i in eachindex(field_names)
        Plots.plot!(plts[i], 2π./ω, getfield(d, field_names[i]), xscale=:log10, 
            yerr = getfield(d_err, field_names[i]),
            yscale= yscale[i], 
            ylabel= "$(units[i])",
            xlabel= "T (s)"; kwargs...
        )
    end
    nothing;
end

function plot_response(plts::Vector{Any}; kwargs...)
    Plots.plot(plts..., layout= (length(plts), 1), size= (400, 400*length(plts)), margin= 3Plots.mm; kwargs...)
end

# for any other survey, we will almost always have two options, either a layer model or model defined at differnt points. Our function does it on a layer model, for a model defined at points, one can simply interpolate and choose small grid spacing for the interpolator, converting the point model into a layer model
"""
    plot_model(m::MTModel; max_depth= m.h[end]*5, kwargs...)
plots a 1D MT model, for now.

## Arguments
* `m` : model to be plotted, contains both the model values and the layer thicknesses

## Keyword arguments
* `max_depth` : max depth to which the last layer should extend, defaults to 5 x the total model thickness. If you want to constrain the plots use `ylim` argument.
* `kwargs` : Controls properties of the plot, splats the argument to `Plots.jl` 
"""
function plot_model(m::MTModel; max_depth= 5*sum(m.h), kwargs...) # add depth lims later?, add kwargs

    h2= cumsum(m.h);

    xlabel, ylabel= get_model_labels(m);
    if max_depth <= h2[end]
        @warn("`max_depth` you provided does not capture the model in entirety. \n Defaulting it to 5 x the total model thickness. \n If you want to constrain the plots use `ylim` argument.") 
        max_depth= 5*h2[end];
    end
    
    Plots.plot([m.m[1], m.m...], [1,h2..., max_depth], scale=:log10, linetype=:steppost,
        yflip= true, xlabel= xlabel, ylabel= ylabel;
        kwargs...
    )
    
end

"""
    plot_model!(plt, m::MTModel; max_depth= 5*cumsum(m.h), kwargs...) 
Plots on the same model with the same functionality as `plot_model`.
"""
function plot_model!(plt, m::MTModel; max_depth= 5*cumsum(m.h), kwargs...) 

    h2= cumsum(m.h);

    xlabel, ylabel= get_model_labels(m);
    if max_depth <= h2[end]
        @warn("`max_depth` you provided does not capture the model in entirety. \n Defaulting it to 2 x the total model thickness. \n If you want to constrain the plots use `ylim` argument.") 
        max_depth= 2*h2[end];
    end
    
    Plots.plot!(plt, [m.m[1], m.m...], [1, h2..., max_depth], scale=:log10, linetype=:steppost,
        yflip= true;
        kwargs...
    )
    
end
