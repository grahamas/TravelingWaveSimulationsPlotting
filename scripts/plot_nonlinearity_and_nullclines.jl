using TravelingWaveSimulations, WilsonCowanModel, TravelingWaveSimulationsPlotting, Simulation73Plotting

include("../_drafts/plot_nonlinearity.jl")

# prototype = get_prototype("full_dynamics_monotonic")
# sim = prototype()

# params = HE2018Params(sim)

prototype = get_prototype("full_dynamics_blocking")
mods = (α=(0.4, 0.7), Aie=0.81, Aei=0.8, firing_θI=0.2, blocking_θI=0.5, save_idxs=nothing, save_on=true, saveat=0.1) 
sim = prototype(; mods...)

@show calculate_fixedpoints.(Ref(sim.model), [0.01, 0.001])
params = get_nullcline_params(sim)

nullcline_sc, nullcline_ly = nullclines(params, 0.01);
display(nullcline_sc)

nonl_scene, nonl_layout = layout_plot(plot_nonlinearity!, sim);
display(nonl_scene)

exec = execute(sim);
exec_sc = exec_heatmap(exec);
display(exec_sc)

for strength in 0.0:0.02:0.1
    @show strength
    sim = prototype(; mods..., stim_strength=strength)
    exec = execute(sim);
    exec_sc = exec_heatmap(exec);
    save(plotsdir("he_heatmap_fp7_str$(strength).png"), exec_sc)
end

