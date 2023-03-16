function plot(m::model; label= false, title= "") # add depth lims later?, add kwargs
    plot([m.ρ[1], m.ρ...], [1, m.h..., m.h[end]*3], linetype=:steppost, scale=:log10, yflip= true, 
        label= label, title= title,
        xlabel= "ρₐ (Ωm)",
        ylabel= "depth (m)"
    )
end
function plot!(m::model; label= false, title= "")
    plot!([m.ρ[1], m.ρ...], [1, m.h..., m.h[end]*3], linetype=:steppost, scale=:log10, yflip= true, 
        label= label, title= title,
        xlabel= "ρₐ (Ωm)",
        ylabel= "depth (m)"
    )
end

function plot(d::response, ω; r::Symbol=:none, label= false, title= "")
    if r===:ρₐ || r===:app_res
        plot(2π./ω, d.ρₐ, scale=:log10,
            label= label, title= title,
            xlabel= "T (s)", ylabel= "ρₐ (Ωm)",
        )
    elseif r===:ϕ || r===:phase
        plot(2π./ω, d.ϕ, xscale=:log10,
            label= label, title= title,
            xlabel= "T (s)", ylabel= "ϕ (⁰)",
        )
    else
        p1= plot(2π./ω, d.ρₐ, scale=:log10,
            label= label, title= title,
            xlabel= "T (s)", ylabel= "ρₐ (Ωm)",
        )
        
        p2= plot(2π./ω, d.ϕ, xscale=:log10,
            label= label, title= title,
            xlabel= "T (s)", ylabel= string(r)* "ϕ (⁰)",
        )
        
        plot!(p1, p2, layout= (2,1), size= (600, 600))
    end
end