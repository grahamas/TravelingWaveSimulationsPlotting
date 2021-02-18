using DrWatson
using TravelingWaveSimulations, WilsonCowanModel, 
    TravelingWaveSimulationsPlotting, 
    Simulation73Plotting, Simulation73
using Dates
using Makie

let example_name = "he_blocking_7fp",
    session_id = "$(Dates.now())",
    example_dir = mkpath(plotsdir("$(example_name)_$(session_id)")),
    stim_strengths = [0.0, 0.02, 0.1],
    figure_resolution = (3000, 1600);
prototype = get_prototype("full_dynamics_blocking")
mods = (α=(0.4, 0.7), 
        Aie=0.81, Aei=0.8, 
        firing_θI=0.2, θI=0.2, blocking_θI=0.5, 
        save_idxs=nothing, save_on=true, saveat=0.1, 
        stop_time=4., x_lattice=256., n_lattice=512) 
sim = prototype(; mods...)

@show calculate_fixedpoints.(Ref(sim.model), [0.01, 0.001])
# params = get_nullcline_params(sim)
# nullcline_scene, nullcline_ly = nullclines(params, 0.01);
# save(joinpath(example_dir, "$(example_name)_nullclines.png"), nullcline_scene)

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

with_theme(simple_theme) do 
    fig = Figure(resolution = figure_resolution)
    fig[1,1] = plot_nonlinearity!(fig, sim; title=nothing)
    fig[2,1] = plot_connectivity!(fig, sim; xlims=[-200, 200], title=nothing)
    fig[1:2,2] = plot_stimulus!(fig, sim)

    execs_by_stim_strength = map(stim_strengths) do stim_strength
        sim = prototype(; mods..., stim_strength=stim_strength)
        execute(sim)
    end
    max_val = maximum(maximum.(execs_by_stim_strength))
    min_val = 0
    color_limits = [min_val, ceil(max_val, sigdigits=2)]

    n_methods_cols = 2
    for (i_exec, exec) in enumerate(execs_by_stim_strength)
        fig[1:2,n_methods_cols+i_exec] = exec_heatmap!(fig, exec; clims=color_limits,
            title = "S = $(stim_strengths[i_exec])", colorbar_width = 0)
    end
    exec_axis = content(fig[1:2,n_methods_cols+1][2,1])
    exec_heatmap = exec_axis.scene.plots |> only
    cbar = fig[2,n_methods_cols+length(stim_strengths)+1] = Colorbar(fig, exec_heatmap, height = exec_axis.height, width=25, vertical=true, ticks=color_limits)#,ticklabels=string.(clims))

    label_a = fig[1, 1, TopLeft()] = Label(fig, "A", textsize = 35,
                #font = noto_sans_bold, 
                halign = :left)
    label_b = fig[2, 1, TopLeft()] = Label(fig, "B", textsize = 35,
                #font = noto_sans_bold, 
                halign = :left)
    label_c = fig[1, 2, TopLeft()] = Label(fig, "C", textsize = 35,
                #font = noto_sans_bold, 
                halign = :center)
    label_d = fig[1, 3, TopLeft()] = Label(fig, "D", textsize = 35,
                #font = noto_sans_bold, 
                halign = :left)

    save(joinpath(example_dir, "$(example_name)_summary.png"), fig) 
end

# with_theme(simple_theme) do 
#     nonlinearity_figure = figure_plot(plot_nonlinearity!, sim; resolution = (800, 800));
#     save(joinpath(example_dir, "$(example_name)_nonlinearity.png"), nonlinearity_figure)

#     connectivity_figure = figure_plot(plot_connectivity!, sim; resolution = (800, 800), xlims=[-200., 200.]);
#     save(joinpath(example_dir, "$(example_name)_connectivity.png"), connectivity_figure)

#     stimulus_figure = figure_plot(plot_stimulus!, sim; resolution = (800, 1200));
#     save(joinpath(example_dir, "$(example_name)_stimulus.png"), stimulus_figure)
# end



# for stim_strength in stim_strengths
#     sim = prototype(; mods..., stim_strength=stim_strength)
#     exec = execute(sim);
#     exec_scene = exec_heatmap(exec);
#     save(joinpath(example_dir, "$(example_name)_heatmap_strength$(stim_strength).png"), exec_scene)
# end

end
