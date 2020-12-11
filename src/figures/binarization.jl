function figure_metasweep_binarization_counts!(scene::Scene, mdb_path::String, 
        (x_sym, y_sym)::Tuple{Symbol,Symbol}, 
        metasweep_syms::NTuple{N,Symbol}, 
        property_sym::Symbol; plots_subdir=nothing, 
        label_translate=default_label_translate,
        errorband_fn=x -> confint(OneSampleTTest(x)) .- mean(x)) where N
    data = TravelingWaveSimulations.load_classifications(mdb_path)[property_sym]

    metasweep_values = map(sym -> axes_keys(data)[dim(data, sym)], metasweep_syms)
    metasweep_value_pairs = product(metasweep_values...)
    getindex_metasweep_dim(arr, idxs) = getindex(arr; Dict(sym => idx for (sym, idx) in zip(metasweep_syms, idxs))...)

    (xs, ys) = axes_keys(_collapse_to_axes(getindex_metasweep_dim(data, first(metasweep_value_pairs)), x_sym, y_sym))

    # FIXME indexing for testing
    binarization_counts = map(metasweep_value_pairs) do metasweep_value_pair
        metasweep_slice = getindex_metasweep_dim(data, metasweep_value_pair)
        total_propagating = count(metasweep_slice)
        total_sims = prod(size(metasweep_slice))
        propagation = (n_propagating=total_propagating, n_total=total_sims)
        flattened_metasweep_slice = _collapse_to_axes(metasweep_slice, x_sym, y_sym)
        binary_segmentation = calc_binary_segmentation(flattened_metasweep_slice)
        if plots_subdir !== nothing
            save_reduce_2d_and_steepest_line_and_histogram((x_sym, y_sym), flattened_metasweep_slice,
                    property_sym, "$metasweep_value_pair", plots_subdir; facet_title="$metasweep_value_pair")
        end
        return binary_segmentation, propagation
    end

    binary_segmentations = [tup[1] for tup in binarization_counts]
    propagations = [tup[2] for tup in binarization_counts]

    set_theme!(LAxis=(textsize=5,), LText=(tellwidth=false, tellheight=false))
    layout = GridLayout()

    layout[1,1] = metasweep_plot!(scene, metasweep_values, [bin_seg.some / sum(bin_seg) for bin_seg in binary_segmentations]; plotted_var_name="sometimes propagating", metasweep_var_names=label_translate.(metasweep_syms))
    layout[1,2] = metasweep_plot!(scene, metasweep_values, [bin_seg.all / sum(bin_seg) for bin_seg in binary_segmentations]; plotted_var_name="always propagating", metasweep_var_names=label_translate.(metasweep_syms))

    layout[2,1] = metasweep_plot!(scene, metasweep_values, [prop.n_propagating / prop.n_total for prop in propagations]; plotted_var_name="overall proportion propagating", metasweep_var_names=label_translate.(metasweep_syms))

    return layout
end