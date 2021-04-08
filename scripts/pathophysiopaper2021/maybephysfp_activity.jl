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
    session_name = "fp_activity",
    session_id = "$(Dates.now())",
    figure_resolution=(700,700)
;

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

put_into(target, source) = target[1:length(source)] = source

for (sym, fp_arr) ∈ [(:monotonic, monotonic_fp_arr), (:blocking, blocking_fp_arr)]
    all_fp_Es = Array{Union{Float64,Missing}}(missing, size(fp_arr)..., 7)
    for mx_idx ∈ Iterators.product(axes(fp_arr)...)
        fps = fp_arr[mx_idx...]
        for fp_idx ∈ 1:length(fps) 
            all_fp_Es[mx_idx..., fp_idx] = first(fps[fp_idx])
        end
    end
    plot_and_save_ax(hist!,
        skipmissing(all_fp_Es) |> collect,
        figure_resolution=figure_resolution,
        plot_name="allfp_$(sym)_excitatory_activity_hist.png",
        title="Histogram of excitatory activity of all FP ($(sym))",
        plots_subdir=plots_subdir
    )
end

# for n_fp in 1:2:7
#     blocking_fp_idx = blocking_fp_count_arr .== n_fp
#     blocking_fp_fps = blocking_fp_arr[blocking_fp_idx]

#     mid_fp = n_fp ÷ 2 + 1
#     blocking_fp_middlefps = map(blocking_fp_fps) do fps
#         fps_sorted_by_E = sort(fps, by=first)
#         fps_sorted_by_E[mid_fp]
#     end

#     with_theme(simple_theme) do

#         plot_and_save_ax(hist!,
#             blocking_fp_middlefps .|> first,
#             figure_resolution=figure_resolution,
#             plot_name="middleblocking_fp_of_$(n_fp)_excitatory_activity_hist.png",
#             title="Histogram of excitatory activity of middle FP (block) (of $n_fp)",
#             plots_subdir=plots_subdir
#         )

#         plot_and_save_ax(hist!,
#             blocking_fp_middlefps .|> x -> x[end],
#             figure_resolution=figure_resolution,
#             plot_name="middlefpblocking_of_db$(n_fp)_inhibitory_activity_hist.png",
#             title="Histogram of inhibitory activity of middle FP (block) (of $n_fp)",
#             plots_subdir=plots_subdir
#         )
#     end # with_theme

#     monotonic_fp_idx = monotonic_fp_count_arr .== n_fp
#     monotonic_fp_fps = monotonic_fp_arr[monotonic_fp_idx]

#     mid_fp = n_fp ÷ 2 + 1
#     monotonic_fp_middlefps = map(monotonic_fp_fps) do fps
#         fps_sorted_by_E = sort(fps, by=first)
#         fps_sorted_by_E[mid_fp]
#     end

#     with_theme(simple_theme) do

#         plot_and_save_ax(hist!,
#             monotonic_fp_middlefps .|> first,
#             figure_resolution=figure_resolution,
#             plot_name="middlemonotonic_fp_of_$(n_fp)_excitatory_activity_hist.png",
#             title="Histogram of excitatory activity of middle FP (mono) (of $n_fp)",
#             plots_subdir=plots_subdir
#         )

#         plot_and_save_ax(hist!,
#             monotonic_fp_middlefps .|> x -> x[end],
#             figure_resolution=figure_resolution,
#             plot_name="middlemonotonic_fp_of_db$(n_fp)_inhibitory_activity_hist.png",
#             title="Histogram of inhibitory activity of middle FP (mono) (of $n_fp)",
#             plots_subdir=plots_subdir
#         )
#     end # with_theme
# end #
end # let

