using TravelingWaveSimulations, WilsonCowanModel, TravelingWaveSimulationsPlotting, Simulation73Plotting,
Simulation73, DrWatson, Dates

#using GLMakie; ext_2d = "png"; GLMakie.activate!()
using CairoMakie; ext_2d = "svg"; CairoMakie.activate!(); 

let example_name = "pathophysio_example",
    session_id = "$(Dates.now())",
    stim_strengths = [0.0, 0.02, 0.1, 0.5, 1.0],
    mods = (
        α=(0.4, 0.7), 
        firing_θI=0.2, blocking_θI=0.5, 
        n_lattice = 512,
        save_idxs=nothing, save_on=true, saveat=0.1
    ),
    figure_resolution=(500,500)
;
example_dir = mkpath(plotsdir("$(example_name)_$(session_id)"))
prototype = get_prototype("full_dynamics_blocking")
sim = prototype(; mods...)

@show length.(calculate_fixedpoints.(Ref(sim.model), [100, 1000]))
params = get_nullcline_params(sim)

plot_and_save(plot_nullclines!, params, 100; 
    plots_subdir=example_dir,
    plot_name="$(example_name)_nullclines.$(ext_2d)",
    figure_resolution=figure_resolution
)

plot_and_save(plot_nonlinearity!, sim; 
    plots_subdir=example_dir,
    plot_name="$(example_name)_nonlinearity.$(ext_2d)",
    figure_resolution=figure_resolution
)

for stim_strength in stim_strengths
    sim = prototype(; mods..., stim_strength=stim_strength)
    exec = execute(sim);
    plot_and_save(exec_heatmap!, exec;
        plots_subdir=example_dir,
        plot_name="$(example_name)_heatmap_strength$(stim_strength).$(ext_2d)"
    ) 
end

end
