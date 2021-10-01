using DrWatson
using TravelingWaveSimulations, WilsonCowanModel, 
    TravelingWaveSimulationsPlotting, 
    Simulation73Plotting, Simulation73
using TravelingWaveSimulationsPlotting: _collapse_to_axes
using Dates
#using GLMakie; ext_2d = "png"; GLMakie.activate!()
using CairoMakie; ext_2d = "svg"; CairoMakie.activate!()
using AxisIndices

# loads blocking_fp_arr and monotonic_fp_arr
sub_A_sweep_lower_bound = 0.5
sub_A_sweep_upper_bound = 1.5
sub_A_range = sub_A_sweep_lower_bound..sub_A_sweep_upper_bound
include(scriptsdir("load/he_sweep_arrs.jl"))



let blocking_fp_arr = blocking_fp_arr[Aee=sub_A_range,Aei=sub_A_range,Aie=sub_A_range,Aii=sub_A_range], 
    monotonic_fp_arr = monotonic_fp_arr[Aee=sub_A_range,Aei=sub_A_range,Aie=sub_A_range,Aii=sub_A_range],
    blocking_fp_count_arr = length.(blocking_fp_arr), 
    monotonic_fp_count_arr = length.(monotonic_fp_arr),
    blocking_prototype_name = "full_dynamics_blocking",
    monotonic_prototype_name = "full_dynamics_monotonic",
    arrows_step=0.05,
    smods = (
        α=(0.4, 0.7), 
        firing_θI=0.2, blocking_θI=0.5, 
        n_lattice = 2,
        save_idxs=nothing, save_on=false, saveat=0.1
    ),
    session_name = "fp_activity_$sub_A_range",
    session_id = "$(Dates.now())",
    figure_resolution=(800,800),
    abbrev_count_label = x -> begin
        if x >= 1000
            try
                "$(Int(x / 1000))K"
            catch
                "$(x / 1000)K"
            end
        else
            "$(Int(x))"
        end
    end,
    bar_theme = Theme(
        fontsize=30,
        Axis = (
            backgroundcolor = :white,
            leftspinevisible = true,
            rightspinevisible = false,
            bottomspinevisible = true,
            topspinevisible = false,
            xgridcolor = :white,
            ygridcolor = :white,
            ytickformat = xs -> abbrev_count_label.(xs)
        )
    ),
    nullcline_theme = Theme(
        fontsize=30,
        Axis = (
            backgroundcolor = :white,
            leftspinevisible = true,
            rightspinevisible = false,
            bottomspinevisible = true,
            topspinevisible = false,
            xgridcolor = :white,
            ygridcolor = :white
        ),
        Lines = (
            linewidth=4.0,
        ),
        Arrows = (
            arrowsize=10, lengthscale=0.017,
            linewidth=2,
            arrowcolor=:black, linecolor=:black,
            colormap=ColorSchemes.Greys_5,
            normalize=false
        ),
        Scatter = (
            markersize=27,
            strokewidth=1
        )
    )  
;

with_theme(bar_theme) do
blocking_stable_fp_arr = filter_stable_fps(blocking_prototype_name, smods, blocking_fp_arr)
monotonic_stable_fp_arr = filter_stable_fps(monotonic_prototype_name, smods, monotonic_fp_arr)


n_sfp_blocking = count_stable_fps(blocking_prototype_name, smods, blocking_fp_arr)
n_sfp_monotonic = count_stable_fps(monotonic_prototype_name, smods, monotonic_fp_arr)

@show count(n_sfp_blocking .== 4)
@show count(n_sfp_monotonic .== 3)
@show count(n_sfp_blocking .== 3 .& n_sfp_monotonic .== 3)

plots_subdir = plotsdir("$(session_name)_$(session_id)")
mkpath(plots_subdir)

blocking_3sfp_si = seizure_index.(blocking_fp_arr[n_sfp_blocking .== 3]) |> collect
@show length(blocking_3sfp_si)
@show count(blocking_3sfp_si .> 0.7) / length(blocking_3sfp_si)
plot_and_save_ax(hist!,
    blocking_3sfp_si,
    figure_resolution=figure_resolution,
    plot_name="blocking_3sfp_si_hist.$(ext_2d)",
    Axis = (
        title="3SFP (blocking)",
        ylabel="count",
        xlabel="max seizure index"
    ),
    plots_subdir=plots_subdir
)

plot_and_save_ax(hist!,
    blocking_fp_arr[n_sfp_blocking .== 3] .|> fps -> maximum(first.(fps)),
    figure_resolution=figure_resolution,
    plot_name="blocking_3sfp_max_excitation_hist.$(ext_2d)",
    Axis = (
        title="3SFP (blocking)",
        ylabel="count",
        xlabel="max E"
    ),
    plots_subdir=plots_subdir
)

monotonic_3sfp_si = seizure_index.(monotonic_fp_arr[n_sfp_monotonic .== 3]) |> collect
@show length(monotonic_3sfp_si)
plot_and_save_ax(hist!,
    monotonic_3sfp_si,
    figure_resolution=figure_resolution,
    plot_name="monotonic_3sfp_si_hist.$(ext_2d)",
    Axis = (
        title="3SFP (monotonic)",
        ylabel="count",
        xlabel="max seizure index"
    ),
    plots_subdir=plots_subdir
)

plot_and_save_ax(hist!,
    monotonic_fp_arr[n_sfp_monotonic .== 3] .|> fps -> maximum(first.(fps)),
    figure_resolution=figure_resolution,
    plot_name="monotonic_3sfp_max_excitation_hist.$(ext_2d)",
    Axis = (
        title="3SFP (monotonic)",
        ylabel="count",
        xlabel="max E"
    ),
    plots_subdir=plots_subdir
)

plot_and_save_ax(hist!,
    seizure_index.(monotonic_fp_arr[n_sfp_blocking .== 3]) |> collect,
    figure_resolution=figure_resolution,
    plot_name="monotonic_former3sfp_si_hist.$(ext_2d)",
    Axis = (
        title="3SFP blocking -> monotonic",
        ylabel="count",
        xlabel="max seizure index"
    ),
    plots_subdir=plots_subdir
)

plot_and_save_ax(hist!,
    monotonic_fp_arr[n_sfp_blocking .== 3] .|> fps -> maximum(first.(fps)),
    figure_resolution=figure_resolution,
    plot_name="monotonic_former3sfp_max_excitation_hist.$(ext_2d)",
    Axis = (
        title="3SFP blocking -> monotonic",
        ylabel="count",
        xlabel="max E"
    ),
    plots_subdir=plots_subdir
)

plot_and_save_ax(hist!,
    n_sfp_monotonic[n_sfp_blocking .== 3],
    figure_resolution=figure_resolution,
    plot_name="monotonic_former3sfp_sfp_count.$(ext_2d)",
    Axis = (
        title="3SFP blocking -> monotonic",
        ylabel="count",
        xlabel="# SFP"
    ),
    plots_subdir=plots_subdir
)

plot_and_save_ax(hist!,
    seizure_index.(monotonic_fp_arr[seizure_index.(blocking_fp_arr) .> 0.71]) |> collect,
    figure_resolution=figure_resolution,
    plot_name="monotonic_formermaxEp_si_hist.$(ext_2d)",
    Axis = (
        title="Maximal EM score blocking -> monotonic",
        ylabel="count",
        xlabel="max seizure index"
    ),
    plots_subdir=plots_subdir
)

plot_and_save_ax(hist!,
    monotonic_fp_arr[seizure_index.(blocking_fp_arr) .> 0.71] .|> fps -> maximum(first.(fps)),
    figure_resolution=figure_resolution,
    plot_name="monotonic_formermaxEp_max_excitation_hist.$(ext_2d)",
    Axis = (
        title="Maximal EM score blocking -> monotonic",
        ylabel="count",
        xlabel="max E"
    ),
    plots_subdir=plots_subdir
)

plot_and_save_ax(hist!,
    monotonic_fp_arr[map(fps -> maximum(first.(fps)), blocking_fp_arr) .> 0.71] .|> fps -> maximum(first.(fps)),
    figure_resolution=figure_resolution,
    plot_name="monotonic_formermaxE_max_excitation_hist.$(ext_2d)",
    Axis = (
        title="Maximal E blocking -> monotonic",
        ylabel="count",
        xlabel="max E"
    ),
    plots_subdir=plots_subdir
)

plot_and_save_ax(hist!,
    monotonic_stable_fp_arr[map(fps -> maximum(first.(fps)), blocking_fp_arr) .> 0.71] .|> fps -> maximum(first.(fps)),
    figure_resolution=figure_resolution,
    plot_name="monotonic_STABLE_formermaxE_max_excitation_hist.$(ext_2d)",
    Axis = (
        title="Maximal E blocking -> monotonic (stable)",
        ylabel="count",
        xlabel="max E"
    ),
    plots_subdir=plots_subdir
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
    @show "$sym all: $(count(.!ismissing.(all_fp_Es)))"
    plot_and_save_ax(hist!,
        skipmissing(all_fp_Es) |> collect,
        figure_resolution=figure_resolution,
        plot_name="allfp_$(sym)_excitatory_activity_hist.$(ext_2d)",
        plots_subdir=plots_subdir,
        Axis = (
            title="All FP ($(sym))",
            ylabel="count",
            xlabel="E"
        )
    )
end

has_found = false
for (sym, fp_arr, prototype_name) ∈ [(:monotonic, monotonic_fp_arr, "full_dynamics_monotonic"), (:blocking, blocking_fp_arr, "full_dynamics_blocking")] #, smods, prototype_name) ∈ [(:monotonic, monotonic_fp_arr, monotonic_nullcline_mods, monotonic_prototype_name), (:blocking, blocking_fp_arr, blocking_nullcline_mods, blocking_prototype_name)]
    prototype = get_prototype(prototype_name)
    stable_fp_Es = Array{Union{Float64,Missing}}(missing, size(fp_arr)..., 7)
    for (mx_idx, (nt, fps)) ∈ zip(CartesianIndices(fp_arr), TravelingWaveSimulationsPlotting.enumerate_nt(fp_arr))
        params = get_nullcline_params(prototype(; smods..., nt...))
        for fp_idx ∈ 1:length(fps)
            fp = fps[fp_idx]
            if TravelingWaveSimulationsPlotting.fixedpoint_is_stable(params, fp)
                stable_fp_Es[mx_idx, fp_idx] = first(fp)
            end
        end
    end
    @show "$sym stable: $(count(.!ismissing.(stable_fp_Es)))"
    plot_and_save_ax(hist!,
        skipmissing(stable_fp_Es) |> collect,
        figure_resolution=figure_resolution,
        plot_name="stablefp_$(sym)_excitatory_activity_hist.$(ext_2d)",
        plots_subdir=plots_subdir,
        Axis = (
            title="Stable FP ($(sym))",
            ylabel="count",
            xlabel="E"
        )
    )
end

for (sym, fp_arr) ∈ [(:monotonic, monotonic_fp_arr), (:blocking, blocking_fp_arr)] 
    all_fp_si = Array{Union{Float64,Missing}}(missing, size(fp_arr)..., 7)
    for mx_idx ∈ Iterators.product(axes(fp_arr)...)
        fps = fp_arr[mx_idx...]
        for fp_idx ∈ 1:length(fps) 
            all_fp_si[mx_idx..., fp_idx] = seizure_index(fps[fp_idx])
        end
    end
    @show "$sym all: $(count(.!ismissing.(all_fp_si)))"
    plot_and_save_ax(hist!,
        skipmissing(all_fp_si) |> collect,
        figure_resolution=figure_resolution,
        plot_name="allfp_$(sym)_seizure_index_hist.$(ext_2d)",
        plots_subdir=plots_subdir,
        Axis = (
            title="All FP ($(sym))",
            ylabel="count",
            xlabel="EI contrast"
        )
    )
end

for (sym, fp_arr, n_sfp_arr, prototype_name) ∈ [(:monotonic, monotonic_fp_arr, n_sfp_monotonic, "full_dynamics_monotonic"), (:blocking, blocking_fp_arr, n_sfp_blocking, "full_dynamics_blocking")] #, smods, prototype_name) ∈ [(:monotonic, monotonic_fp_arr, monotonic_nullcline_mods, monotonic_prototype_name), (:blocking, blocking_fp_arr, blocking_nullcline_mods, blocking_prototype_name)]
    prototype = get_prototype(prototype_name)
    stable_fp_Es = Array{Union{Float64,Missing}}(missing, size(fp_arr)..., 7)
    for (mx_idx, (nt, fps)) ∈ zip(CartesianIndices(fp_arr), TravelingWaveSimulationsPlotting.enumerate_nt(fp_arr))
        params = get_nullcline_params(prototype(; smods..., nt...))
        for fp_idx ∈ 1:length(fps)
            fp = fps[fp_idx]
            if TravelingWaveSimulationsPlotting.fixedpoint_is_stable(params, fp)
                stable_fp_Es[mx_idx, fp_idx] = seizure_index(fp)
            end
        end
    end
    @show "$sym stable: $(count(.!ismissing.(stable_fp_Es)))"
    plot_and_save_ax(hist!,
        skipmissing(stable_fp_Es) |> collect,
        figure_resolution=figure_resolution,
        plot_name="stablefp_$(sym)_seizure_index_hist.$(ext_2d)",
        plots_subdir=plots_subdir,
        Axis = (
            title="Stable FP ($(sym))",
            ylabel="count",
            xlabel="EI contrast"
        )
    )
    
    with_theme(nullcline_theme) do
    example_3sfp_idxs = findall(n_sfp_arr .== 3)
    example_3sfp_idx = example_3sfp_idxs[end ÷ 2]
    example_3sfp_coord = TravelingWaveSimulationsPlotting.get_coordinate(fp_arr, example_3sfp_idx)
    @show "$sym 3SFP: $example_3sfp_coord"
    plot_and_save(plot_nullclines!,
        get_nullcline_params(prototype(; smods..., example_3sfp_coord...));
        plot_name="nullclines_3sfp_$(sym).$(ext_2d)",
        arrows_step=arrows_step,
        plots_subdir=plots_subdir
    )
    if sym == :blocking
        example_4sfp_idxs = findall(n_sfp_arr .== 4)
        example_4sfp_idx = example_4sfp_idxs[end ÷ 4]
        example_4sfp_coord = TravelingWaveSimulationsPlotting.get_coordinate(fp_arr, example_4sfp_idx)
        @show "4FP: $example_4sfp_coord"
        plot_and_save(plot_nullclines!,
            get_nullcline_params(prototype(; smods..., example_4sfp_coord...));
            plot_name="nullclines_4sfp_$(sym).$(ext_2d)",
            arrows_step=arrows_step,
            plots_subdir=plots_subdir
        )

        example_7fp_idxs = findall(blocking_fp_count_arr[:] .== 7)
        example_7fp_idx = example_7fp_idxs[end]
        example_7fp_coord = TravelingWaveSimulationsPlotting.get_coordinate(fp_arr, example_7fp_idx)
        @show "7FP: $example_7fp_coord"
        plot_and_save(plot_nullclines!,
            get_nullcline_params(prototype(; smods..., example_7fp_coord...));
            plot_name="nullclines_7fp_$(sym).$(ext_2d)",
            arrows_step=arrows_step,
            plots_subdir=plots_subdir
        )
    end
    example_both_3sfp_idxs = findall(n_sfp_monotonic .== 3 .& n_sfp_blocking .== 3)
    example_both_3sfp_idx = example_both_3sfp_idxs[end]
    example_both_3sfp_coord = TravelingWaveSimulationsPlotting.get_coordinate(fp_arr, example_both_3sfp_idx)
    @show "Both: $example_both_3sfp_coord"
    plot_and_save(plot_nullclines!,
        get_nullcline_params(prototype(; smods..., example_both_3sfp_coord...));
        plot_name="nullclines_both_3sfp_$(sym).$(ext_2d)",
        arrows_step=arrows_step,
        plots_subdir=plots_subdir
    )
end
end

for (sym, fp_arr, prototype_name) ∈ [(:blocking, blocking_fp_arr, "full_dynamics_blocking")] 
    prototype = get_prototype(prototype_name)
    mid_7fp_stable_Es = Array{Union{Float64,Missing}}(missing, size(fp_arr)..., 7)
    mid_7fp_stable_si = Array{Union{Float64,Missing}}(missing, size(fp_arr)..., 7)
    mid_7fp_Es = Array{Union{Float64,Missing}}(missing, size(fp_arr)..., 7)
    mid_7fp_si = Array{Union{Float64,Missing}}(missing, size(fp_arr)..., 7)
    for (mx_idx, (nt, fps)) ∈ zip(CartesianIndices(fp_arr), TravelingWaveSimulationsPlotting.enumerate_nt(fp_arr))
        params = get_nullcline_params(prototype(; smods..., nt...))
        if length(fps) == 7
            fp = fps[4]
            if TravelingWaveSimulationsPlotting.fixedpoint_is_stable(params, fp)
                mid_7fp_stable_Es[mx_idx, 4] = first(fp)
                mid_7fp_stable_si[mx_idx, 4] = seizure_index(fp)
            end
            mid_7fp_Es[mx_idx, 4] = first(fp)
            mid_7fp_si[mx_idx, 4] = seizure_index(fp)
        end
    end
    @show "$sym stable: $(count(.!ismissing.(mid_7fp_Es)))"
    @show "$sym stable: $(count(.!ismissing.(mid_7fp_stable_Es)))"
    plot_and_save_ax(hist!,
        skipmissing(mid_7fp_si) |> collect,
        figure_resolution=figure_resolution,
        plot_name="mid7fp_$(sym)_seizure_index_hist.$(ext_2d)",
        Axis = (
            title="mid-7FP ($sym)",
            ylabel="count",
            xlabel="EI contrast"
        ),
        plots_subdir=plots_subdir
    )
    plot_and_save_ax(hist!,
        skipmissing(mid_7fp_stable_si) |> collect,
        figure_resolution=figure_resolution,
        plot_name="mid7fp_stable_$(sym)_seizure_index_hist.$(ext_2d)",
        Axis = (
            title="stable mid-7FP ($sym)",
            ylabel="count",
            xlabel="EI contrast"
        ),
        plots_subdir=plots_subdir
    )
     plot_and_save_ax(hist!,
        skipmissing(mid_7fp_Es) |> collect,
        figure_resolution=figure_resolution,
        plot_name="mid7fp_$(sym)_excitatory_hist.$(ext_2d)",
        Axis = (
            title="mid-7FP ($sym)",
            ylabel="count",
            xlabel="E"
        ),
        plots_subdir=plots_subdir
    )
    plot_and_save_ax(hist!,
        skipmissing(mid_7fp_stable_Es) |> collect,
        figure_resolution=figure_resolution,
        plot_name="mid7fp_stable_$(sym)_excitatory_hist.$(ext_2d)",
        Axis = (
            title="stable mid-7FP ($sym)",
            ylabel="count",
            xlabel="E"
        ),
        plots_subdir=plots_subdir
    )
end

end # with_theme
end # let

