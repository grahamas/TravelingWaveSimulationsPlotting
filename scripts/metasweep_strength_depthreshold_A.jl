using TravelingWaveSimulationsPlotting, TravelingWaveSimulations

metasweep_A_fpath = TravelingWaveSimulations.get_recent_simulation_data_path(joinpath(homedir(), "data", "ring_normed_blocking", "strength_depthreshold_A"))

# Definitely propagating
plot_and_save(figure_metasweep_binarization_counts!, metasweep_A_fpath;
    to_binarize_axes_syms=(:Aei, :Aee), 
    independent_axes_syms=(:blocking_θI, :stim_strength), 
    property_sym=:has_propagation,
    scene_resolution=(2400,2400), 
    sweep_name="strength_depthreshold_A")


                           