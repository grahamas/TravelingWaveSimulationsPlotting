using TravelingWaveSimulations, WilsonCowanModel, TravelingWaveSimulationsPlotting, Simulation73Plotting,
        DrWatson

include("../_drafts/plot_nonlinearity.jl")

# prototype = get_prototype("full_dynamics_monotonic")
# sim = prototype()

# params = HE2018Params(sim)
using Dates

let example_name = "he_2018",
    session_id = "$(Dates.now())",
    example_dir = mkpath(plotsdir("$(example_name)_$(session_id)")),
    stim_strengths = [0.0, 0.02, 0.1, 0.5, 1.0];
prototype = get_prototype("harris_ermentrout")
mods = (#α=(0.4, 0.7), Aie=0.81, Aei=0.8, firing_θI=0.2, blocking_θI=0.5,
        stim_duration=20., stop_time=20., τ=(1.0,0.1),
        Sii=0.001, Sei=0.001, See=1.0, Sie=1.0, x_lattice=512.,
        save_idxs=nothing, save_on=true, saveat=0.1) 
sim = prototype(; mods...)

@show calculate_fixedpoints.(Ref(sim.model), [0.01, 0.001])
params = get_nullcline_params(sim)


nullcline_scene, nullcline_ly = nullclines(params, 0.01);
save(joinpath(example_dir, "$(example_name)_nullclines.png"), nullcline_scene)

nonl_scene, nonl_layout = layout_plot(plot_nonlinearity!, sim);
save(joinpath(example_dir, "$(example_name)_nonlinearity.png"), nonl_scene)

for stim_strength in stim_strengths
    sim = prototype(; mods..., stim_strength=stim_strength)
    exec = execute(sim);
    # @show exec.solution.u .|> u -> population(u, 1) |> pop -> pop[begin:min(length(pop), 500)] |> maximum
    # @show exec.solution.u .|> u -> population(u, 2) |> pop -> pop[begin:min(length(pop), 500)] |> maximum
    exec_scene = exec_heatmap(exec);
    save(joinpath(example_dir, "$(example_name)_heatmap_strength$(stim_strength).png"), exec_scene)
end

end
