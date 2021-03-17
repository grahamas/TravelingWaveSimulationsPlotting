using DrWatson
using TravelingWaveSimulations, WilsonCowanModel, 
    TravelingWaveSimulationsPlotting, 
    Simulation73Plotting, Simulation73
using TravelingWaveSimulationsPlotting: _collapse_to_axes
using Dates
using Makie
using AxisIndices

# loads blocking_fp_arr and monotonic_fp_arr
include(scriptsdir("load/he_sweep_arrs.jl"))


let blocking_fp_arr = blocking_fp_arr, 
    monotonic_fp_arr = monotonic_fp_arr,
    blocking_fp_count_arr = length.(blocking_fp_arr), 
    monotonic_fp_count_arr = length.(monotonic_fp_arr),
    session_name = "updown_vs_epileptic",
    session_id = "$(Dates.now())",
    colorbar_width = 15,
    figure_resolution=(700,700)
;

possible_fp_counts = 0:1:7 # FIXME should be odd only
possible_visible_axes = [(:Aee, :Aei), (:Aie, :Aii), (:Aee, :Aie), (:Aei, :Aii), (:Aee, :Aii), (:Aei, :Aie)]

plots_subdir = plotsdir("$(session_name)_$(session_id)")
mkpath(plots_subdir)

simple_theme = Theme(
    linewidth = 20.0,
    fontsize=24,
    Axis = (
        backgroundcolor = :white,
        leftspinevisible = true,
        rightspinevisible = false,
        bottomspinevisible = true,
        topspinevisible = false,
        xgridcolor = :white,
        ygridcolor = :white,
    )
)

with_theme(simple_theme) do
    mono_5fp = monotonic_fp_count_arr .== 5
    block_7fp = blocking_fp_count_arr .== 7
    both_5_and_7 = mono_5fp .& block_7fp
    prop_5_to_7 = count(both_5_and_7) / count(mono_5fp) * 100

    shouldbezeros_not_5_but_7 = block_7fp .& .!mono_5fp
    shouldbezero_count = count(shouldbezeros_not_5_but_7)
    @info "sanity: 7s not from 5s: $(shouldbezero_count)"
    if shouldbezero_count != 0
        plot_and_save(axisarray_heatmap!,
            _collapse_to_axes(shouldbezeros_not_5_but_7, :Aee, :Aei);
            figure_resolution=figure_resolution,
            colorbar_width=colorbar_width,
            plot_name="SANITY_CHECK.png",
            title="Blocking-7FP that do NOT come from mono-5FP",
            plots_subdir=plots_subdir
        )
    end
    @info "proportion 5FP-mono becoming 7FP-block: $(prop_5_to_7)%"

    plot_and_save(axisarray_heatmap!,
        _collapse_to_axes(both_5_and_7, :Aee, :Aei);
        figure_resolution=figure_resolution,
        colorbar_width=colorbar_width,
        plot_name="mono5fp_block7fp_overlap.png",
        title="Overlap between mono-5FP and blocking-7FP",
        plots_subdir=plots_subdir
    )

    monotonic_epilepsy_metric_arr = epilepsy_metric.(monotonic_fp_arr)
    blocking_epilepsy_metric_arr = epilepsy_metric.(blocking_fp_arr)

    plot_and_save_ax(hist!,
        monotonic_epilepsy_metric_arr[mono_5fp];
        figure_resolution=figure_resolution,
        colorbar_width=colorbar_width,
        plot_name="mono5fp_epilepsy_metric_hist.png",
        title="Histogram of epilepsy metric for mono-5FP",
        plots_subdir=plots_subdir
    )

    plot_and_save_ax(hist!,
        blocking_epilepsy_metric_arr[block_7fp];
        figure_resolution=figure_resolution,
        colorbar_width=colorbar_width,
        plot_name="block7fp_epilepsy_metric_hist.png",
        title="Histogram of epilepsy metric for block-7FP",
        plots_subdir=plots_subdir
    )

end # with_theme
end # let