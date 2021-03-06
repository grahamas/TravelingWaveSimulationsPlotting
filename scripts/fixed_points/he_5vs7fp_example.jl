using DrWatson
using TravelingWaveSimulations, WilsonCowanModel, 
    TravelingWaveSimulationsPlotting, 
    Simulation73Plotting, Simulation73
using TravelingWaveSimulationsPlotting: _collapse_to_axes, calculate_fixedpoints
using Dates
using Makie

include(projectdir("_drafts/grouped_bar_plot.jl"))

# loads blocking_fp_arr and monotonic_fp_arr
include(scriptsdir("load/he_sweep_arrs.jl"))

let example_name = "he_5vs7fp",
    session_id = "$(Dates.now())",
    example_dir = mkpath(plotsdir("$(example_name)_$(session_id)")),
    stim_strengths = [0.0, 0.02, 0.1],
    figure_resolution = (3000, 1600),
    colorbar_width = 25;

blocking_fp_count_arr=length.(blocking_fp_arr)
monotonic_fp_count_arr=length.(monotonic_fp_arr)

mods = (α=(0.4, 0.7), 
    Aie=0.81, Aei=0.8, 
    firing_θI=0.2, θI=0.2, blocking_θI=0.5, 
    save_idxs=nothing, save_on=true, saveat=0.1,
    x_lattice=256., n_lattice=512,
    stim_radius = 1.0,
    stop_time = 4
) 
simple_theme = Theme(
    linewidth = 20.0,
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
nonl_type_strings = ["monotonic", "blocking"]
#nonl_type_strings = ["blocking"]

with_theme(simple_theme) do 
    fig = Figure(resolution = figure_resolution)
    for (i_nonl_type, nonl_type_string) ∈ enumerate(nonl_type_strings)
        prototype = get_prototype("full_dynamics_$(nonl_type_string)");
        sim = prototype(; mods...)
        params = get_nullcline_params(sim)

        # test that the FP counter isn't unstable near this model
        test_axis_lengths = [50, 73, 100, 853]
        test_fp_counts = length.(calculate_fixedpoints.(Ref(sim.model), test_axis_lengths))
        @show test_fp_counts
        a_test_is_different = !all(y->y==test_fp_counts[1], test_fp_counts)
        @show a_test_is_different
        
        if a_test_is_different || !isodd(test_fp_counts[begin]) 
            @show calculate_fixedpoints.(Ref(sim.model), test_axis_lengths) .|> length
            fig[1,1] = plot_nullclines!(fig, params, test_axis_lengths[1+0])
            fig[2,1] = plot_nullclines!(fig, params, test_axis_lengths[1+1])
            fig[1,2] = plot_nullclines!(fig, params, test_axis_lengths[1+2])
            fig[2,2] = plot_nullclines!(fig, params, test_axis_lengths[1+3])
            display(fig)
            error("assert failed.")
        end
        fig[i_nonl_type, 1] = plot_nullclines!(fig, params, test_axis_lengths[1])
        execs_by_stim_strength = [
            execute(prototype(; mods..., stim_strength=stim_strength)) 
            for stim_strength ∈ stim_strengths
        ]
        fig[i_nonl_type, 1 .+ (1:length(stim_strengths))] = 
            isolimit_exec_heatmaps!(
                fig, execs_by_stim_strength
                ; titles = ["S = $stim_str" for stim_str ∈ stim_strengths]
        );
    end # for nonl_type

    length(nonl_type_strings) != 2 && @warn "Figure construction incorrectly assumes exactly two nonlinearity types."
    summary_fig_layout = GridLayout()

    # Histogram of #fixed points
    fp_hist_ax = let fig = fig
        n_fps = 0:7
        blocking_fp_count_hist = [log10(count(blocking_fp_count_arr .== x)) for x in n_fps]
        monotonic_fp_count_hist = [log10(count(monotonic_fp_count_arr .== x)) for x in n_fps]
        max_log = vcat(blocking_fp_count_hist, monotonic_fp_count_hist) |> maximum |> mx -> ceil(Int, mx)

        ax = Makie.Axis(fig)
        gbp = groupedbarplot!(ax, n_fps, [monotonic_fp_count_hist, blocking_fp_count_hist])
        gbp_labels = ["mono", "block"]
        tightlimits!(ax)
        xlims!(ax, n_fps[begin]-0.5,n_fps[end]+0.5)
        ylims!(ax, 0, max_log)
        ax.xticks=n_fps
        ax.yticks=0:max_log
        ax.ytickformat = xs -> [x > 0 ? "10^$(Int(x))" : "$x" for x in xs]
        ax.xlabel = "# fixed points"
        ax.ylabel = "# models"
        hidespines!(ax, :t, :r)
        hidedecorations!(ax, label=false, ticklabels=false, ticks=false)
        axislegend(ax, gbp, gbp_labels)
        ax
    end # let fp_hist_ax

    fp_diff_arr = blocking_fp_count_arr .- monotonic_fp_count_arr

    # Heatmap of difference #FP(block - mono) == target_fp_diff
    fp_diff_count_ax = let fig = fig,
            target_fp_diff = 2,
            visible_diff_axes = (:Aee, :Aei),
            fp_diff_arr = fp_diff_arr,
            colorbar_width = colorbar_width;
        smushed_diff_count_arr = _collapse_to_axes(fp_diff_arr .== target_fp_diff, visible_diff_axes...)
        axisarray_heatmap!(fig, smushed_diff_count_arr;
            colorbar_width = colorbar_width,
            colorbar_label = "∝FP(block - mono) == 2"
        )
    end # let fp_diff_count_ax
    hidedecorations!(content(fp_diff_count_ax[1,1]))

    # Heatmap of #FP(block) == target_fp
    fp_count_ax = let fig = fig,
            target_fp = 7,
            visible_diff_axes = (:Aee, :Aei),
            fp_count_arr = blocking_fp_count_arr,
            colorbar_width = colorbar_width;
        smushed_count_arr = _collapse_to_axes(fp_count_arr .== target_fp, visible_diff_axes...)
        axisarray_heatmap!(fig, smushed_count_arr;
            colorbar_width=colorbar_width,
            colorbar_label="∝FP(block) == 7"
        )
    end # let fp_diff_count_ax

    summary_fig_layout[:v] = [fp_diff_count_ax,
        fp_count_ax,
        fp_hist_ax]

    fig[1:2, 2+length(stim_strengths)] = summary_fig_layout

    label_a = fig[1, 1, TopLeft()] = Label(fig, "A", textsize = 35,
                #font = noto_sans_bold, 
                halign = :left)
    label_b = fig[1, 2, TopLeft()] = Label(fig, "B", textsize = 35,
                #font = noto_sans_bold, 
                halign = :left)
    label_c = fig[2, 1, TopLeft()] = Label(fig, "C", textsize = 35,
                #font = noto_sans_bold, 
                halign = :left)
    label_d = fig[2, 2, TopLeft()] = Label(fig, "D", textsize = 35,
                #font = noto_sans_bold, 
                halign = :left)
    label_e = fig[1:2, 2+length(stim_strengths)][1,1, TopLeft()] = 
        Label(fig, "E", textsize = 35,
                #font = noto_sans_bold, 
                halign = :left)
    label_f = fig[1:2, 2+length(stim_strengths)][3,1, TopLeft()] = 
        Label(fig, "F", textsize = 35,
                #font = noto_sans_bold, 
                halign = :left)

    save(
        joinpath(example_dir, 
            "$(example_name).png"
        ), 
        fig
    )
end # with_theme
end # let
