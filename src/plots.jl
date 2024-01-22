struct utils
    ρₐ::Symbol
    ϕ::Symbol
end

# most geophysical surveys generally have the same units across, what would be the point of passing units then?

function prepare_plot(d::response, ω::Vector{T}; field_names::Vector{Symbol}= [k for k ∈ fieldnames(typeof(d))], kwargs...) where T <: Union{Float32, Float64} 
    
    units= utils(:Ωm, :ᵒ);
    yscale= utils(:log10, :identity)
    plts= [];

    for i in eachindex(field_names)
        pi= Plots.plot(2π./ω, getfield(d, field_names[i]), xscale=:log10, 
            yscale= getfield(yscale, field_names[i]),
            ylabel= "$(field_names[i]) ($(getfield(units, field_names[i])))",
            xlabel= "T (s)"; kwargs...
        )
        push!(plts, pi);
    end

    return plts
end

function prepare_plot!(plts::Vector{Any}, d::response, ω; field_names::Vector{Symbol}= [k for k ∈ fieldnames(typeof(d))], kwargs...)

    units= utils(:Ωm, :ᵒ);
    yscale= utils(:log10, :identity)

    for i in eachindex(field_names)
        Plots.plot!(plts[i], 2π./ω, getfield(d, field_names[i]), xscale=:log10, 
            yscale= getfield(yscale, field_names[i]),
            ylabel= "$(field_names[i]) ($(getfield(units, field_names[i])))",
            xlabel= "T (s)"; kwargs...
        )
    end
    nothing;
end

function plot_response(plts::Vector{Any}; kwargs...)
    Plots.plot(plts..., layout= (length(plts), 1), size= (400, 400*length(plts)), margin= 3Plots.mm; kwargs...)
end


function get_model_labels(m::model) # can be defined for different types of models of other surveys
    return "ρ (Ωm)", "Depth (m)"; 
end

# for any other survey, we will almost always have two options, either a layer model or model defined at differnt points. Our function does it on a layer model, for a model defined at points, one can simply interpolate and choose small grid spacing for the interpolator, converting the point model into a layer model

function plot_model(m::model; max_depth= m.h[end]*5, kwargs...) # add depth lims later?, add kwargs

    h2= cumsum(m.h);

    xlabel, ylabel= get_model_labels(m);
    if max_depth <= h2[end]
        @warn("`max_depth` you provided does not capture the model in entirety. \n Defaulting it to 2 x the total model thickness. \n If you want to constrain the plots use `ylim` argument.") 
        max_depth= 2*h2[end];
    end
    
    Plots.plot([m.m[1], m.m...], [1,h2..., max_depth], scale=:log10, linetype=:steppost,
        yflip= true, xlabel= xlabel, ylabel= ylabel;
        kwargs...
    )
    
end

function plot_model!(plt, m::model; max_depth= m.h[end]*5, kwargs...) 

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
