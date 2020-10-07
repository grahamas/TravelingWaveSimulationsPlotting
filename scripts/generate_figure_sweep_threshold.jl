using TravelingWaveSimulationsPlotting, TravelingWaveSimulations

sweep_threshold_A_fpath = TravelingWaveSimulations.get_recent_simulation_data_path(joinpath(homedir(), "data", "ring_blocking", "depblock_A"))

# Definitely propagating
save_sweep_threshold(sweep_threshold_A_fpath, (:Aei, :Aee), :blocking_θI, 
                            :has_propagation)
                           