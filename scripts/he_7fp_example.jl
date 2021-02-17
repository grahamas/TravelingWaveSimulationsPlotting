using TravelingWaveSimulations, WilsonCowanModel, TravelingWaveSimulationsPlotting, Simulation73Plotting,
        DrWatson

# prototype = get_prototype("full_dynamics_monotonic")
# sim = prototype()

# params = HE2018Params(sim)
using Dates

let example_name = "he_blocking_7fp",
    session_id = "$(Dates.now())",
    example_dir = mkpath(plotsdir("$(example_name)_$(session_id)")),
    stim_strengths = [0.0, 0.02, 0.1, 0.5, 1.0];
prototype = get_prototype("full_dynamics_blocking")
mods = (α=(0.4, 0.7), Aie=0.81, Aei=0.8, firing_θI=0.2, blocking_θI=0.5, save_idxs=nothing, save_on=true, saveat=0.1) 
sim = prototype(; mods...)

@show calculate_fixedpoints.(Ref(sim.model), [0.01, 0.001])
params = get_nullcline_params(sim)

# nullcline_scene, nullcline_ly = nullclines(params, 0.01);
# save(joinpath(example_dir, "$(example_name)_nullclines.png"), nullcline_scene)
using Makie
ggplot_theme = Theme(
    linewidth = 200,
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

with_theme(ggplot_theme) do 
    nonlinearity_figure = figure_plot(plot_nonlinearity!, sim; resolution = (800, 800));
    save(joinpath(example_dir, "$(example_name)_nonlinearity.png"), nonlinearity_figure)

    connectivity_figure = figure_plot(plot_connectivity!, sim; resolution = (800, 800), xlims=[-200., 200.]);
    save(joinpath(example_dir, "$(example_name)_connectivity.png"), connectivity_figure)
end



# for stim_strength in stim_strengths
#     sim = prototype(; mods..., stim_strength=stim_strength)
#     exec = execute(sim);
#     exec_scene = exec_heatmap(exec);
#     save(joinpath(example_dir, "$(example_name)_heatmap_strength$(stim_strength).png"), exec_scene)
# end

end
