using TravelingWaveSimulationsPlotting, TravelingWaveSimulations, StaticArrays

# Generate figure of a connectivity sweep with examples

threshold_sweep_monotonic_A_fpath = get_recent_simulation_data_path(joinpath(homedir(), "data", "ring_blocking", "depblock_A"), -1)
@show threshold_sweep_monotonic_A_fpath


# Definitely propagating
save_threshold_sweep((:Aei, :Aee), (:Aie, :Aii), [
                          (Aee=200.0, Aei=50.0, Aii=200, Aie=50.0),
                          (Aee=150.0, Aei=150.0, Aii=50, Aie=115.0),
                          (Aee=50.0, Aei=200.0, Aii=200, Aie=50.0)
                         ],  
                            fcmb_monotonic_A_fpath, fcmb_blocking_A_fpath, 
                            :has_propagation, "contrast_monotonic_blocking_A")
                             

save_figure_example_contrast_monotonic_blocking_all((:Sei, :See), (:Sie, :Sii), [
                          (See=14.0, Sei=110.0, Sii=14., Sie=110.0),
                          (See=14.0, Sei=14.0, Sii=110., Sie=14.0),
                          #(See=30.0, Sei=50.0, Sii=50., Sie=90.0),
                          (See=50.0, Sei=40.0, Sii=40., Sie=50.0)
                         ],  
                            fcmb_monotonic_S_fpath, fcmb_blocking_S_fpath, 
                            :has_propagation, "contrast_monotonic_blocking_S")
                           