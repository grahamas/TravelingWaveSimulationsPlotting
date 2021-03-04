
using DrWatson
using Dates
using TravelingWaveSimulationsPlotting: _collapse_to_axes
using TravelingWaveSimulationsPlotting
using StaticArrays, AxisIndices, NamedDims

if !@isdefined(refresh_sweep_arrs)
    refresh_sweep_arrs = false
end

if !@isdefined(blocking_fp_arr) || !@isdefined(monotonic_fp_arr) || refresh_sweep_arrs
    A_range = 0.1:0.1:1.5
    sweeping_mods = (Aee=A_range, Aei=A_range, Aie=A_range, Aii=A_range)
    static_mods = (
        α=(0.4, 0.7), 
        firing_θI=0.2, blocking_θI=0.5, 
        n_lattice = 2,
        save_idxs=nothing, save_on=true, saveat=0.1
    ) 
    file, filename = produce_or_load(datadir(), []; 
        prefix = "blocking_fp_arr",
        force = refresh_sweep_arrs
    ) do c
        blocking_fp_arr = sweep_calculate_fixedpoints(
            "full_dynamics_blocking", 
            static_mods,
            sweeping_mods,
            ; 
            dx = 0.01
        )
        return @dict(blocking_fp_arr)
    end
    @unpack blocking_fp_arr = file
    file, filename = produce_or_load(datadir(), []; 
        prefix = "monotonic_fp_arr",
        force = refresh_sweep_arrs
    ) do c
        monotonic_fp_arr = sweep_calculate_fixedpoints(
            "full_dynamics_monotonic", 
            static_mods,
            sweeping_mods,
            ; 
            dx = 0.01
        )
        return @dict(monotonic_fp_arr)
    end
    @unpack monotonic_fp_arr = file
end


let fp_arrs = [("blocking", blocking_fp_arr), 
                ("monotonic", monotonic_fp_arr)],
    session_name = "he_fp_sweeps",
    session_id = "$(Dates.now())",
    colorbar_width = 5

possible_fp_counts = 0:1:7 # FIXME should be odd only
possible_visible_axes = [(:Aee, :Aei), (:Aie, :Aii), (:Aee, :Aie), (:Aei, :Aii), (:Aee, :Aii), (:Aei, :Aie)]

plots_subdir = plotsdir("$(session_name)_$(session_id)")
mkpath(plots_subdir)

for (nonl_type, fp_arr) in fp_arrs
    fp_count_arr = length.(fp_arr)
    fp_epilepsy_arr = epilepsy_metric.(fp_arr)
    for visible_axes in possible_visible_axes
        smushed_count_arr = _collapse_to_axes(fp_count_arr, visible_axes...)
        plot_and_save(axisarray_heatmap!, smushed_count_arr; 
            colorbar_width = colorbar_width,
            plot_name="fpcount_$(nonl_type)_$(join(visible_axes, "_")).png",
            plots_subdir=plots_subdir, colorrange=(0,7))
        smushed_epilepsy_arr = _collapse_to_axes(fp_epilepsy_arr, visible_axes...)
        plot_and_save(axisarray_heatmap!, smushed_count_arr; 
            colorbar_width = colorbar_width,
            plot_name="fpepilepsy_$(nonl_type)_$(join(visible_axes, "_")).png",
            plots_subdir=plots_subdir, colorrange=(-1,1))
        for fp_count in possible_fp_counts
            smushed_count_arr = _collapse_to_axes(fp_count_arr .== fp_count, 
                                            visible_axes...)
            plot_and_save(axisarray_heatmap!, smushed_count_arr;
                colorbar_width = colorbar_width,
                plot_name="fpcountequals$(fp_count)_$(nonl_type)_$(join(visible_axes, "_")).png",
                plots_subdir=plots_subdir)
        end
    end
end

nonl_types = [tup[1] for tup in fp_arrs]
blocking_idx = findfirst(nonl_types .== "blocking")
monotonic_idx = findfirst(nonl_types .== "monotonic")
fp_diff_arr = fp_count_arrs[blocking_idx][2] .- fp_count_arrs[monotonic_idx][2]

for visible_axes in possible_visible_axes
    smushed_arr = _collapse_to_axes(fp_diff_arr, visible_axes...)
    plot_and_save(axisarray_heatmap!, smushed_arr; 
            colorbar_width = colorbar_width,
            plot_name="fpdiff_$(join(visible_axes, "_")).png",
            plots_subdir=plots_subdir)
    for fp_diff in minimum(fp_diff_arr):1:maximum(fp_diff_arr)
        smushed_count_arr = _collapse_to_axes(fp_diff_arr .== fp_diff, 
                                        visible_axes...)
        plot_and_save(axisarray_heatmap!, smushed_count_arr;
            colorbar_width = colorbar_width,
            plot_name="fpdiffequals$(fp_diff)_$(join(visible_axes, "_")).png",
            plots_subdir=plots_subdir)
    end
end
end #let

