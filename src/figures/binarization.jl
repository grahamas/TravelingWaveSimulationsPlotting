export figure_binarization_counts!, multiple_figures_binarization_counts

function multiple_figures_binarization_counts(mdb_path::String;
        multiplot_dimension::Symbol,
        property_sym::Symbol,
        scene_resolution=(1600,1600),
        kwargs...)
    mdb = TravelingWaveSimulations.load_classifications(mdb_path);
    data = mdb[property_sym];

    multiplot_values = keys.(axes(data))[findfirst(AxisIndices.NamedDims.names(data) .== multiplot_dimension)]

    map(multiplot_values) do multiplot_value
        scene, layout = layoutscene(resolution=scene_resolution)
        layout[1,1] = figure_binarization_counts!(scene, getindex(data; Dict(multiplot_dimension => multiplot_value)...); kwargs...)
        (scene, layout, "$(string(multiplot_dimension))=$(multiplot_value)")
    end
end

function multiple_separate_subfigures_binarization_counts(mdb_path::String;
        multiplot_dimension::Symbol,
        property_sym::Symbol,
        scene_resolution=(1600,1600),
        kwargs...)
    mdb = TravelingWaveSimulations.load_classifications(mdb_path);
    data = mdb[property_sym];

    multiplot_values = keys.(axes(data))[findfirst(AxisIndices.NamedDims.names(data) .== multiplot_dimension)]

    n_binarizations = 3
    binarizations_types = ["some", "always", "overall"]
    sep_subfigs = map(multiplot_values) do multiplot_value
        scene_layouts = [layoutscene(resolution=scene_resolution) for _ in 1:n_binarizations]
        figure_binarization_counts!(scene_layouts, getindex(data; Dict(multiplot_dimension => multiplot_value)...); kwargs...)
        [(scene_layout..., "$(string(multiplot_dimension))=$(multiplot_value)_$(binarizations_types[i])") for (i, scene_layout) in enumerate(scene_layouts)]
    end
    return vcat(sep_subfigs...)
end
export multiple_separate_subfigures_binarization_counts

function figure_binarization_counts!(target::Union{Scene,Vector{Tuple}}, mdb_path::String; 
        property_sym::Symbol,
        subslice=(;),
        plots_subdir=nothing, 
        kwargs...)
    """
        Load database from mdb_path, grabbing AxisArray stored at `property_sym` as `data`

        Take `subslice` of `data`, if specified

        For every coordinate along `independent_axes` specified by `independent_axes_syms`, Construct N_binarize-dim array of proportion of simulations that propagated

        Then binarize and count those proportions, producing N_independent-dim plots showing the "binarized" counts as a function of `independent_axes_syms`

        TODO: binarize is a bad name. Here I actually trinarize into `none`, `some`, and `all` propagating for each bin
    """
    mdb = TravelingWaveSimulations.load_classifications(mdb_path);
    data = mdb[property_sym];
    data = getindex(data; subslice...);

    figure_binarization_counts!(scene, data; kwargs...)
end

function figure_binarization_counts!(scene::Scene, data::AbstractArray; 
        to_binarize_axes_syms::NTuple{N_binarize, Symbol}, 
        independent_axes_syms::NTuple{N_independent, Symbol}, 
        plots_subdir=nothing, 
        label_translate=default_label_translate,
        errorband_fn=x -> confint(OneSampleTTest(x)) .- mean(x),
        session_name=nothing, session_id=nothing) where {N_binarize, N_independent}
    """
        Load database from mdb_path, grabbing AxisArray stored at `property_sym` as `data`

        Take `subslice` of `data`, if specified

        For every coordinate along `independent_axes` specified by `independent_axes_syms`, Construct N_binarize-dim array of proportion of simulations that propagated

        Then binarize and count those proportions, producing N_independent-dim plots showing the "binarized" counts as a function of `independent_axes_syms`

        TODO: binarize is a bad name. Here I actually trinarize into `none`, `some`, and `all` propagating for each bin
    """
    
    independent_axes_values = map(sym -> axes_keys(data)[findfirst(AxisIndices.NamedDims.names(data) .== sym)], independent_axes_syms)
    independent_axes_grid = product(independent_axes_values...)
    getindex_on_independent_axes(arr, idxs) = getindex(arr; Dict(sym => idx for (sym, idx) in zip(independent_axes_syms, idxs))...)

    binarization_counts_by_independent_axes = map(independent_axes_grid) do independent_coord
        slice_to_binarize = getindex_on_independent_axes(data, independent_coord)
        total_propagating = count(slice_to_binarize)
        total_sims = prod(size(slice_to_binarize))
        propagation = (n_propagating=total_propagating, n_total=total_sims)
        flattened_slice_to_binarize = _collapse_to_axes(slice_to_binarize, to_binarize_axes_syms...)
        binary_segmentation = TravelingWaveSimulationsPlotting.calc_binary_segmentation(flattened_slice_to_binarize)
        return binary_segmentation, propagation
    end

    binary_segmentations = [tup[1] for tup in binarization_counts_by_independent_axes]
    propagations = [tup[2] for tup in binarization_counts_by_independent_axes]

    set_theme!(Axis=(textsize=5,), Label=(tellwidth=false, tellheight=false))
    layout = GridLayout()

    layout[1,1] = metasweep_plot!(scene, independent_axes_values, [bin_seg.some / sum(bin_seg) for bin_seg in binary_segmentations]; 
            plotted_var_name="sometimes propagating", 
            lims = (0.,1.), 
            metasweep_var_names=label_translate.(independent_axes_syms))
    layout[1,2] = metasweep_plot!(scene, independent_axes_values, [bin_seg.all / sum(bin_seg) for bin_seg in binary_segmentations]; plotted_var_name="always propagating", 
            lims = (0.,1.), 
            metasweep_var_names=label_translate.(independent_axes_syms))

    layout[2,1] = metasweep_plot!(scene, independent_axes_values, [prop.n_propagating / prop.n_total for prop in propagations]; plotted_var_name="overall proportion propagating", 
            lims = (0.,1.), 
            metasweep_var_names=label_translate.(independent_axes_syms))

    return layout
end


function figure_binarization_counts!(scene_layouts::AbstractVector{<:Tuple}, data::AbstractArray; 
        to_binarize_axes_syms::NTuple{N_binarize, Symbol}, 
        independent_axes_syms::NTuple{N_independent, Symbol}, 
        plots_subdir=nothing, 
        label_translate=default_label_translate,
        errorband_fn=x -> confint(OneSampleTTest(x)) .- mean(x),
        session_name=nothing, session_id=nothing) where {N_binarize, N_independent}
    """
        Load database from mdb_path, grabbing AxisArray stored at `property_sym` as `data`

        Take `subslice` of `data`, if specified

        For every coordinate along `independent_axes` specified by `independent_axes_syms`, Construct N_binarize-dim array of proportion of simulations that propagated

        Then binarize and count those proportions, producing N_independent-dim plots showing the "binarized" counts as a function of `independent_axes_syms`

        TODO: binarize is a bad name. Here I actually trinarize into `none`, `some`, and `all` propagating for each bin
    """
    
    independent_axes_values = map(sym -> axes_keys(data)[findfirst(AxisIndices.NamedDims.names(data) .== sym)], independent_axes_syms)
    independent_axes_grid = product(independent_axes_values...)
    getindex_on_independent_axes(arr, idxs) = getindex(arr; Dict(sym => idx for (sym, idx) in zip(independent_axes_syms, idxs))...)

    binarization_counts_by_independent_axes = map(independent_axes_grid) do independent_coord
        slice_to_binarize = getindex_on_independent_axes(data, independent_coord)
        total_propagating = count(slice_to_binarize)
        total_sims = prod(size(slice_to_binarize))
        propagation = (n_propagating=total_propagating, n_total=total_sims)
        flattened_slice_to_binarize = _collapse_to_axes(slice_to_binarize, to_binarize_axes_syms...)
        binary_segmentation = TravelingWaveSimulationsPlotting.calc_binary_segmentation(flattened_slice_to_binarize)
        return binary_segmentation, propagation
    end

    binary_segmentations = [tup[1] for tup in binarization_counts_by_independent_axes]
    propagations = [tup[2] for tup in binarization_counts_by_independent_axes]

    # set_theme!(Axis=(textsize=5,), Label=(tellwidth=false, tellheight=false))

    bin_seg_total = binary_segmentations[1].total
    prop_total = propagations[1].n_total
    @assert all([bin_seg.total == bin_seg_total for bin_seg in binary_segmentations])
    @assert all([prop.n_total == prop_total for prop in propagations])
    scene_layouts[1][2][1,1] = metasweep_plot!(scene_layouts[1][1], independent_axes_values, [float(bin_seg.some) for bin_seg in binary_segmentations]; 
            plotted_var_name="some propagating (of $bin_seg_total)", 
            lims = (0.,bin_seg_total), 
            metasweep_var_names=label_translate.(independent_axes_syms))
    scene_layouts[2][2][1,1] = metasweep_plot!(scene_layouts[2][1], independent_axes_values, [float(bin_seg.all) for bin_seg in binary_segmentations]; plotted_var_name="all propagating (of $bin_seg_total)", 
            lims = (0.,bin_seg_total), 
            metasweep_var_names=label_translate.(independent_axes_syms))

    scene_layouts[3][2][1,1] = metasweep_plot!(scene_layouts[3][1], independent_axes_values, [float(prop.n_propagating) for prop in propagations]; plotted_var_name="count propagating (of $prop_total)", 
            lims = (0., prop_total), 
            metasweep_var_names=label_translate.(independent_axes_syms))

    return scene_layouts
end