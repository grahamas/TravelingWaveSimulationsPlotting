using Makie, Colors

@recipe(GroupedBarPlot, x, groups) do scene
    default_theme(scene, BarPlot)
end

Makie.conversion_trait(::Type{GroupedBarPlot}) = Makie.NoConversion()

function Makie.plot!(p::GroupedBarPlot)
    widths = lift(p.width, p[1], p[2]) do width, x, groups
        n_x, n_groups = length(x), length(groups)
        extra_space = sum(diff(x)) / n_x # half of this on either end
        group_span = (x[end] - x[begin] + extra_space) / n_x
        if width === Makie.automatic
            bars_span = group_span * 0.8
            bar_span = bars_span / n_groups
            [bar_span for group in groups]
        elseif width isa Number
            @assert width > 0
            @assert width * n_groups <= group_span
            [width for group in groups]
        elseif width isa AbstractVector
            @assert length(width) == length(groups)
            @assert sum(width) <= group_span
            width
        else
            error("Unsupported type for GroupedBarPlot attribute `width`")
        end
    end

    group_xs = lift(p[1], p[2], widths) do x, groups, widths
        bars_span = sum(widths)
        # first bar starts at: x - (bars_span/2) + (bar_span/2)
        group_xs = map(enumerate(groups)) do (i_group, group)
            # first bar start + preceding bars + half current bar
            midpoint_offset = -bars_span / 2 + sum(widths[1:i_group-1]) + widths[i_group]/2
            x .+ midpoint_offset
        end
        return group_xs
    end

    colors = lift(p.color, p[2], p.parent.attributes[:backgroundcolor], p.attributes[:strokecolor]) do color, groups, bgcolor, fgcolor
        n_groups = length(groups)
        if color isa AbstractArray && length(color) == n_groups
            color
        else
            distinguishable_colors(n_groups, parse.(Colorant, [bgcolor, fgcolor]), dropseed=true)
        end
    end

    p = lift(group_xs, p[2], colors, widths) do group_xs, group_ys, colors, widths
        for (group_x, group_y, color, width) ∈ zip(group_xs, group_ys, colors, widths)
            barplot!(p, group_x, group_y; p.attributes..., color=color, width=width)
        end
        p
    end

    p
end

Makie.Legend(fig_or_scene, grouped_bar::GroupedBarPlot, args...; kwargs...) = Legend(fig_or_scene, grouped_bar.plots, args...; kwargs...)
