using TravelingWaveSimulationsPlotting, TravelingWaveSimulations
using Dates

session_name = "wider_strength_depthreshold_A"

metasweep_A_fpath = TravelingWaveSimulations.get_recent_simulation_data_path(joinpath(homedir(), "data", "ring_normed_blocking", session_name))

# Definitely propagating
session_id = "$(Dates.now())"

plot_and_save_many(multiple_averaged_to_axes,
    metasweep_A_fpath;
    axes=[(:Aei, :Aee), (:Aie, :Aii)],
    multi_param_syms=[:blocking_θI, :stim_strength],
    property_sym=:has_propagation,
    scene_resolution=(1200, 400),
    session_name=session_name,
    session_id=session_id,
)
                           
# plot_and_save_many(multiple_separate_subfigures_binarization_counts, metasweep_A_fpath;
#     to_binarize_axes_syms=(:Aei, :Aee), 
#     independent_axes_syms=(:blocking_θI,),
#     multiplot_dimension=:stim_strength,
#     property_sym=:has_propagation,
#     scene_resolution=(400,400), 
#     session_name=session_name,
#     session_id=session_id)
                           
# plot_and_save_many(multiple_separate_subfigures_binarization_counts, metasweep_A_fpath;
#     to_binarize_axes_syms=(:Aei, :Aee), 
#     independent_axes_syms=(:stim_strength,),
#     multiplot_dimension=:blocking_θI,
#     property_sym=:has_propagation,
#     scene_resolution=(400,400), 
#     session_name=session_name,
#     session_id=session_id)