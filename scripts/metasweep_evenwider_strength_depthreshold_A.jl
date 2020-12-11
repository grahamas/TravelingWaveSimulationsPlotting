using TravelingWaveSimulationsPlotting, TravelingWaveSimulations

sweep_name = "evenwider_strength_depthreshold_A"

metasweep_A_fpath = TravelingWaveSimulations.get_recent_simulation_data_path(joinpath(homedir(), "data", "ring_blocking", sweep_name))

# plot_and_save(figure_metasweep_binarization_counts!, metasweep_A_fpath, (:Aei, :Aee), (:blocking_θI, :stim_strength), 
#                             :has_propagation; scene_resolution=(2400,2400), 
#                             sweep_name=sweep_name)
                           
# plot_and_save(figure_example!, [(metasweep_A_fpath, (Aee=600., Aei=160., Aie=300., Aii=300., blocking_θI=12., stim_strength=13.)),
#             (metasweep_A_fpath, (Aee=300., Aei=160., Aie=300., Aii=300., blocking_θI=12., stim_strength=13.)),
#             (metasweep_A_fpath, (Aee=150., Aei=160., Aie=300., Aii=300., blocking_θI=12., stim_strength=13.)),
#             (metasweep_A_fpath, (Aee=0., Aei=160., Aie=300., Aii=300., blocking_θI=12., stim_strength=13.))
#             ]; sweep_name="examples_stim13_θ12");
        
# plot_and_save(figure_example!, [(metasweep_A_fpath, (Aee=600., Aei=160., Aie=300., Aii=300., blocking_θI=12., stim_strength=2.)),
#             (metasweep_A_fpath, (Aee=300., Aei=160., Aie=300., Aii=300., blocking_θI=12., stim_strength=2.)),
#             (metasweep_A_fpath, (Aee=150., Aei=160., Aie=300., Aii=300., blocking_θI=12., stim_strength=2.)),
#             (metasweep_A_fpath, (Aee=0., Aei=160., Aie=300., Aii=300., blocking_θI=12., stim_strength=2.))
#             ]; sweep_name="examples_stim2_θ12");

plot_and_save(figure_example!, [(metasweep_A_fpath, (Aee=150., Aei=120., Aie=300., Aii=300., blocking_θI=12., stim_strength=str)) for str in 6.6:0.02:6.8]; sweep_name="examples_Aee150_θ12_6.6-6.8stimsweep", animate=false);

plot_and_save(figure_example!, [(metasweep_A_fpath, (Aee=200., Aei=120., Aie=150., Aii=150., blocking_θI=10., stim_strength=str)) for str in 2.:15.]; sweep_name="examples_Aee200_θ10_2-15stimsweep", animate=false);

plot_and_save(figure_example!, [(metasweep_A_fpath, (Aee=95., Aei=120., Aie=150., Aii=150., blocking_θI=10., stim_strength=str)) for str in 2.:15.]; sweep_name="examples_Aee95_θ10_2-15stimsweep", animate=false);

"done"