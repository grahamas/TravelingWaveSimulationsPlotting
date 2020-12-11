using TravelingWaveSimulationsPlotting, TravelingWaveSimulations

let sweep_name = "reduced_strength_depthreshold_A"

    metasweep_A_fpath = TravelingWaveSimulations.get_recent_simulation_data_path(joinpath(homedir(), "data", "ring_normed_blocking", sweep_name))

    # Definitely propagating
    plot_and_save(figure_metasweep_binarization_counts!, metasweep_A_fpath, (:Aei, :Aee), (:blocking_θI, :stim_strength), 
                                :has_propagation; scene_resolution=(2400,2400), 
                                sweep_name=sweep_name)

end


                           