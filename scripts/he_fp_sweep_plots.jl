
using DrWatson
using Dates
using TravelingWaveSimulationsPlotting: _collapse_to_axes

if !@isdefined(blocking_fp_arr)
    @warn "Calculating blocking_fp_arr..."
    include(scriptsdir("he_blocking_fp_sweep.jl"))
    @warn "done."
end

if !@isdefined(mono_fp_arr)
    @warn "Calculating monotonic_fp_arr..."
    include(scriptsdir("he_monotonic_fp_sweep.jl"))
    @warn "done."
end

let fp_arrs = [("blocking", blocking_fp_arr), 
                ("monotonic", mono_fp_arr)],
    session_name = "he_fp_sweeps",
    session_id = "$(Dates.now())",
    colorbar_width = 5

possible_fp_counts = 1:2:7
possible_visible_axes = [(:Aee, :Aei), (:Aie, :Aii), (:Aee, :Aie), (:Aei, :Aii), (:Aee, :Aii), (:Aei, :Aie)]

plots_subdir = plotsdir("$(session_name)_$(session_id)")
mkpath(plots_subdir)

for (nonl_type, fp_arr) in fp_arrs
    for visible_axes in possible_visible_axes
        smushed_arr = _collapse_to_axes(fp_arr, visible_axes...)
        plot_and_save(axisarray_heatmap!, smushed_arr, colorbar_width;
            plot_name="fpcount_$(nonl_type)_$(join(visible_axes, "_")).png",
            plots_subdir=plots_subdir, colorrange=(0,7))
        for fp_count in possible_fp_counts
            smushed_count_arr = _collapse_to_axes(fp_arr .== fp_count, 
                                            visible_axes...)
            plot_and_save(axisarray_heatmap!, smushed_count_arr, colorbar_width;
                plot_name="fpcountequals$(fp_count)_$(nonl_type)_$(join(visible_axes, "_")).png",
                plots_subdir=plots_subdir)
        end
    end
end

nonl_types = [tup[1] for tup in fp_arrs]
blocking_idx = findfirst(nonl_types .== "blocking")
monotonic_idx = findfirst(nonl_types .== "monotonic")
fp_diff_arr = fp_arrs[blocking_idx][2] .- fp_arrs[monotonic_idx][2]

for visible_axes in possible_visible_axes
    smushed_arr = _collapse_to_axes(fp_diff_arr, visible_axes...)
    plot_and_save(axisarray_heatmap!, smushed_arr, colorbar_width;
            plot_name="fpdiff_$(join(visible_axes, "_")).png",
            plots_subdir=plots_subdir)
    for fp_diff in minimum(fp_diff_arr):1:maximum(fp_diff_arr)
        smushed_count_arr = _collapse_to_axes(fp_diff_arr .== fp_diff, 
                                        visible_axes...)
        plot_and_save(axisarray_heatmap!, smushed_count_arr, colorbar_width;
            plot_name="fpdiffequals$(fp_diff)_$(join(visible_axes, "_")).png",
            plots_subdir=plots_subdir)
    end
end
end #let

