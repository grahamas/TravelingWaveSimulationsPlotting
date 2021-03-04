# # Depolarization block meta-sweep
# 

using TravelingWaveSimulationsPlotting

metasweep_A_fpath = joinpath(homedir(), "data", "ring_blocking", "depblock_A", "2020-09-01T13:13:14.552_v1.0-282-g4f19b1f")

# Definitely propagating
save_metasweep_sigmoid_fit(metasweep_A_fpath, (:Aei, :Aee), :blocking_θI, 
                            :has_propagation)
                           