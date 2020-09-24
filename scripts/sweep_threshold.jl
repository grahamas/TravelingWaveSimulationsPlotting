# # Depolarization block meta-sweep
# 

using TravelingWaveSimulationsPlotting

sweep_threshold_A_fpath = joinpath(homedir(), "data", "ring_blocking", "depblock_A", "2020-09-01T13:13:14.552_v1.0-282-g4f19b1f")

# Definitely propagating
save_sweep_threshold(sweep_threshold_A_fpath, (:Aei, :Aee), :blocking_θI, 
                            :has_propagation)
                           