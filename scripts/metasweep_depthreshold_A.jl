using TravelingWaveSimulationsPlotting, TravelingWaveSimulations
using Statistics: std, mean
using HypothesisTests: OneSampleTTest, confint


metasweep_A_fpath = TravelingWaveSimulations.get_recent_simulation_data_path(joinpath(homedir(), "data", "ring_blocking", "depthreshold_A"))

# Definitely propagating
save_metasweep_sigmoid_fit(metasweep_A_fpath, (:Aei, :Aee), (:blocking_θI,), 
                            :has_propagation; scene_resolution=(2400, 2400),
                            errorband_fn = x -> confint(OneSampleTTest(x)) .- mean(x) |> first |> abs)
                           