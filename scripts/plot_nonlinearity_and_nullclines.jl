using TravelingWaveSimulations, WilsonCowanModel, TravelingWaveSimulationsPlotting, Simulation73Plotting

include("../_drafts/plot_nonlinearity.jl")

# prototype = get_prototype("full_dynamics_monotonic")
# sim = prototype()

# params = HE2018Params(sim)

prototype = get_prototype("full_dynamics_blocking")
sim = prototype(; α=(1.08, 1.15), Aie=0.81, Aei=0.8, firing_θI=0.37, blocking_θI=0.65, save_idxs=nothing, save_on=true, saveat=0.1)

params = wcm_nullcline_params(sim)

nullcline_sc, nullcline_ly = nullclines(params, 0.01);

# nonl_scene, nonl_layout = layout_plot(plot_nonlinearity!, sim);

# exec = execute(sim)
# exec_heatmap(exec)