function axisarray_heatmap!(scene::Scene, data::AbstractAxisArray, ax_labels=Union{Tuple,Nothing}, 
        colorbar_width::Union{Nothing,Int}=nothing; hide_y=false)
    sweep_ax = LAxis(scene)
    x, y = axes_keys(data)
    heatmap = heatmap!(sweep_ax, x, y, get_data(data), colorrange=(0,1))
    #tightlimits!(sweep_ax)
    if !(ax_labels isa Nothing)
        (x_sym, y_sym) = ax_labels
        sweep_ax.xlabel = string(x_sym)
        if hide_y
            hideydecorations!(sweep_ax)
        else
            sweep_ax.ylabel = string(y_sym)
        end
    end

    sublayout = GridLayout()
    sublayout[1,1] = sweep_ax
    if colorbar_width !== nothing && abs(-(extrema(data)...)) != 0
        sublayout[1,2] = LColorbar(scene, heatmap, width=colorbar_width)
    end

    return sublayout
end

function axisarray_heatmap!(scene::Scene, data::NamedAxisArray, args...; kwargs...)
    axisarray_heatmap!(scene, data.data, _namedaxisarray_names(data), args...; kwargs...)
end



