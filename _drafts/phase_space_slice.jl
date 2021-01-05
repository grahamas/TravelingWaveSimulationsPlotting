using IterTools, TravelingWaveSimulations, TravelingWaveSimulationsPlotting, AxisIndices

using TravelingWaveSimulationsPlotting: _collapse_to_axes

mdb_path = joinpath(homedir(), "data", "ring_normed_blocking", "wider_strength_depthreshold_A")
to_binarize_axes_syms = (:Aee, :Aei)
independent_axes_syms = (:stim_strength,)#, :blocking_θI)
subslice = (blocking_θI = 9.0,)


mdb = TravelingWaveSimulations.load_classifications(mdb_path |> get_recent_simulation_data_path);
data = mdb[:has_propagation];
data = getindex(data; subslice...);

independent_axes_values = map(sym -> axes_keys(data)[findfirst(AxisIndices.NamedDims.names(data) .== sym)], independent_axes_syms)
independent_grid = product(independent_axes_values...)
getindex_on_independent_axes(arr, idxs) = getindex(arr; Dict(sym => idx for (sym, idx) in zip(independent_axes_syms, idxs))...)

binarization_counts_by_independent_axes = map(independent_grid) do independent_coord
    slice_to_binarize = getindex_on_independent_axes(data, independent_coord)
    total_propagating = count(slice_to_binarize)
    total_sims = prod(size(slice_to_binarize))
    propagation = (n_propagating=total_propagating, n_total=total_sims)
    flattened_slice_to_binarize = _collapse_to_axes(slice_to_binarize, to_binarize_axes_syms...)
    binary_segmentation = TravelingWaveSimulationsPlotting.calc_binary_segmentation(flattened_slice_to_binarize)
    return binary_segmentation, propagation
end