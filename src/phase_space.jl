

function plot_binary_segmentation!(scene::Scene, data)
    segmented = calc_binary_segmentation(data)
    segmented_ax = MakieLayout.Axis(scene)
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
                                                   kwargs...)
    whole_ensemble_data = TravelingWaveSimulations.load_ExecutionClassifications(AbstractArray, fpath)[property_sym]
    data = _collapse_to_axes(whole_ensemble_data, x_sym, y_sym)
    reduce_2d_and_steepest_line_and_histogram!(scene, data, property_sym; kwargs...)
end



function reduce_2d_and_steepest_line_and_histogram!(
                                scene::Scene, 
                                data::AbstractArray,
                                property_sym::Symbol; 
                                facet_title, titlesize=20, hide_y=false,
                                colorbar_width=nothing)
    layout = GridLayout()

    title_facet = Label(scene, facet_title, textsize=titlesize, tellwidth=false)
    sweep_sublayout = axisarray_heatmap!(scene, data; colorbar_width = colorbar_width)
    reduction_sublayout = plot_reduction!(scene, slice, data, get_content(sweep_sublayout, 1, 1))
    if isnothing(reduction_sublayout)
        reduction_sublayout = MakieLayout.Axis(scene)
    end
    segmented_ax = plot_binary_segmentation!(scene, data)

    layout[:v] = [title_facet, sweep_sublayout, reduction_sublayout, segmented_ax]

    return layout
end

function save_reduce_2d_and_steepest_line_and_histogram((x_sym, y_sym),
        fpath::String,
        property_sym, 
        suffix="", root_path=plotsdir(); scene_resolution=(300, 900), kwargs...)
    scene, layout = layoutscene(resolution=scene_resolution)
    layout[1,1] = reduce_2d_and_steepest_line_and_histogram!(scene,
            (x_sym, y_sym),
            fpath,
            property_sym; kwargs...)
    fname = "reduce_2d_and_steepest_line_and_histogram_$(x_sym)_$(y_sym)$(length(suffix) > 0 ? "_$(suffix)" : "").png"
    mkpath(root_path)
    save_path = joinpath(root_path, fname)
    @info "saving $(save_path)"
    Makie.save(save_path, scene)
end

function save_reduce_2d_and_steepest_line_and_histogram((x_sym, y_sym),
        data::AbstractArray,
        property_sym, 
        suffix="", root_path=plotsdir(); scene_resolution=(300, 900), kwargs...)
    scene, layout = layoutscene(resolution=scene_resolution)
    layout[1,1] = reduce_2d_and_steepest_line_and_histogram!(scene,
            data,
            property_sym; kwargs...)
    fname = "reduce_2d_and_steepest_line_and_histogram_$(x_sym)_$(y_sym)$(length(suffix) > 0 ? "_$(suffix)" : "").png"
    mkpath(root_path)
    save_path = joinpath(root_path, fname)
    @info "saving $(save_path)"
    Makie.save(save_path, scene)
end