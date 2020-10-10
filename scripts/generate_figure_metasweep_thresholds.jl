using TravelingWaveSimulationsPlotting, TravelingWaveSimulations

metasweep_A_fpath = TravelingWaveSimulations.get_recent_simulation_data_path(joinpath(homedir(), "data", "ring_blocking", "dep_threshold_A"))

# Definitely propagating
save_metasweep(metasweep_A_fpath, (:Aei, :Aee), :blocking_θI, 
                            :has_propagation)
                           