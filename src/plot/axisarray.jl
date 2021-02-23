_round_extrema((x1, x2); sigdigits) = (floor(x1, sigdigits=sigdigits), ceil(x2, sigdigits=sigdigits))
function axisarray_heatmap!(fig::Figure, data::AxisArray, 
        ax_labels=Union{Tuple,Nothing}, 
        ; colorbar_width::Union{Nothing,Int}=nothing,
        colorrange=_round_extrema(extrema(get_data(data)), sigdigits=2),
        hide_y=false,
        colorbar_label=""
    )
    sweep_ax = AbstractPlotting.Axis(fig)
    x, y = axes_keys(data)
    heatmap = heatmap!(sweep_ax, x, y, get_data(data), colorrange=colorrange)
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
        sublayout[1,2] = Colorbar(fig, heatmap, 
        width=colorbar_width, label=colorbar_label,
        ticks=[colorrange...])
    end

    return sublayout
end

function axisarray_heatmap!(fig::Figure, data::NamedAxisArray, args...; kwargs...)
    axisarray_heatmap!(fig, data.data, _namedaxisarray_names(data), args...; kwargs...)
end



