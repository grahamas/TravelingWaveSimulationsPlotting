

function plot_binary_segmentation!(scene::Scene, data)
    segmented = calc_binary_segmentation(data)
    segmented_ax = LAxis(scene)
    segment_names = string.([keys(segmented)...])
    barplot!(segmented_ax, [values(segmented)...])
    ylims!(segmented_ax, 0, 1)
    segmented_ax.xticks =( 1:3, segment_names)
    segmented_ax.xticklabelrotation = 0.0
    return segmented_ax
end

function get_content(layout, rows, cols, side = MakieLayout.GridLayoutBase.Inner())

    span = MakieLayout.GridLayoutBase.Span(
        MakieLayout.GridLayoutBase.to_ranges(layout, rows, cols)...
    )

    only(filter(layout.content) do c
        c.span == span && c.side == side
    end).content
end

function reduce_2d_and_steepest_line_and_histogram!(scene::Scene, 
                                                   (x_sym, y_sym)::Tuple{Symbol,Symbol},
                                                   fpath::String,
                                                   property_sym::Symbol; 
                                                   facet_title, titlesize=20, hide_y=false,
                                                   colorbar_width=nothing)
    prototype_name, sim_params = read_params_from_data_path(fpath)
    whole_ensemble_data = TravelingWaveSimulations.load_ExecutionClassifications(AbstractArray, fpath)[property_sym]

    data = _collapse_to_axes(whole_ensemble_data, x_sym, y_sym)

    layout = GridLayout()

    title_facet = LText(scene, facet_title, textsize=titlesize, tellwidth=false)
    sweep_sublayout = axisarray_heatmap!(scene, data, colorbar_width)
    reduction_sublayout = plot_max_gradient!(scene, slice, data, get_content(sweep_sublayout, 1, 1))
    segmented_ax = plot_binary_segmentation!(scene, data)

    layout[:v] = [title_facet, sweep_sublayout, reduction_sublayout, segmented_ax]

    return layout
end


function save_reduce_2d_and_steepest_line_and_histogram!((x_sym, y_sym),
        property_sym, 
        unique_id=""; kwargs...)
    scene, layout = layout, scene()
    layout[1,1] = reduce_2d_and_steepest_line_and_histogram!(scene,
            (x_sym, y_sym),
            property_sym; kwargs...)
    fname = "reduce_2d_and_steepest_line_and_histogram_$(x_sym)_$(y_sym).png"
    mkpath(plotsdir(unique_id))
    @info "saving $(plotsdir(unique_id,fname))"
    Makie.save(plotsdir(unique_id,fname), scene)
end