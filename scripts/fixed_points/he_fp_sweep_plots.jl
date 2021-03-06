
using DrWatson
using Dates
using TravelingWaveSimulationsPlotting: _collapse_to_axes
using TravelingWaveSimulationsPlotting
using StaticArrays, AxisIndices, NamedDims
using Makie 

# loads blocking_fp_arr and monotonic_fp_arr
include(scriptsdir("load/he_sweep_arrs.jl"))

let fp_arrs = [("blocking", blocking_fp_arr), 
                ("monotonic", monotonic_fp_arr)],
    session_name = "he_fp_sweeps",
    session_id = "$(Dates.now())",
    colorbar_width = 15,
    figure_resolution=(700,700)

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
for (nonl_type, fp_arr) in fp_arrs
    fp_count_arr = length.(fp_arr)
    fp_epilepsy_arr = seizure_index.(fp_arr)
    for visible_axes in possible_visible_axes
        smushed_count_arr = _collapse_to_axes(fp_count_arr, visible_axes...)
        plot_and_save(axisarray_heatmap!, smushed_count_arr;  
            figure_resolution=figure_resolution,
            colorbar_width = colorbar_width,
            plot_name="fpcount_$(nonl_type)_$(join(visible_axes, "_")).png",
            plots_subdir=plots_subdir, colorrange=(0,7))
        smushed_epilepsy_arr = _collapse_to_axes(fp_epilepsy_arr, visible_axes...)
        plot_and_save(axisarray_heatmap!, smushed_count_arr;  
            figure_resolution=figure_resolution,
            colorbar_width = colorbar_width,
            plot_name="fpepilepsy_$(nonl_type)_$(join(visible_axes, "_")).png",
            plots_subdir=plots_subdir, colorrange=(-1,1))
        for fp_count in possible_fp_counts
            smushed_count_arr = _collapse_to_axes(fp_count_arr .== fp_count, 
                                            visible_axes...)
            plot_and_save(axisarray_heatmap!, smushed_count_arr; 
                figure_resolution=figure_resolution,
                colorbar_width = colorbar_width,
                plot_name="fpcountequals$(fp_count)_$(nonl_type)_$(join(visible_axes, "_")).png", title="$(nonl_type) models with #FP = $(fp_count)",
                plots_subdir=plots_subdir)
        end
    end
end

nonl_types = [tup[1] for tup in fp_arrs]
blocking_idx = findfirst(nonl_types .== "blocking")
monotonic_idx = findfirst(nonl_types .== "monotonic")
fp_diff_arr = length.(fp_arrs[blocking_idx][2]) .- length.(fp_arrs[monotonic_idx][2])

for visible_axes in possible_visible_axes
    smushed_arr = _collapse_to_axes(fp_diff_arr, visible_axes...)
    plot_and_save(axisarray_heatmap!, smushed_arr;  
            figure_resolution=figure_resolution,
            colorbar_width = colorbar_width,
            plot_name="fpdiff_$(join(visible_axes, "_")).png",
            plots_subdir=plots_subdir)
    for fp_diff in minimum(fp_diff_arr):1:maximum(fp_diff_arr)
        smushed_count_arr = _collapse_to_axes(fp_diff_arr .== fp_diff, 
                                        visible_axes...)
        plot_and_save(axisarray_heatmap!, smushed_count_arr; 
            figure_resolution=figure_resolution,
            colorbar_width = colorbar_width,
            plot_name="fpdiffequals$(fp_diff)_$(join(visible_axes, "_")).png",
            plots_subdir=plots_subdir,
            title="#FP(block-mono) = $(fp_diff)")
    end
end
end #with_theme
end #let

