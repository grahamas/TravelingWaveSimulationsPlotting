using TravelingWaveSimulationsPlotting, TravelingWaveSimulations

metasweep_A_fpath = "/home/graham/data/ring_blocking/depblock_A/2020-10-06T22:03:12.176_v1.0-304-g11c6074_dirty"

# Definitely propagating
save_metasweep_sigmoid_fit(metasweep_A_fpath, (:Aei, :Aee), :blocking_θI, 
                            :has_propagation)
                           