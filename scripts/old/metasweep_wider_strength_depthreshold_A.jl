using TravelingWaveSimulationsPlotting, TravelingWaveSimulations
using Dates

session_name = "wider_strength_depthreshold_A"

metasweep_A_fpath = TravelingWaveSimulations.get_recent_simulation_data_path(joinpath(homedir(), "data", "ring_normed_blocking", session_name))

# Definitely propagating
session_id = "$(Dates.now())"

plot_and_save_many(multiple_figure_examples,
    metasweep_A_fpath,
    [
        (Aee=120., Aei=200., Aie=125., Aii=125., blocking_θI=9., stim_strength=2.), # dep-block non-prop. (sus.)
        (Aee=120., Aei=200., Aie=125., Aii=125., blocking_θI=9., stim_strength=3.), # dep-block mid prop. (sus.)
        (Aee=120., Aei=200., Aie=125., Aii=125., blocking_θI=9., stim_strength=20.), # dep-block prop. (sus.)
        (Aee=120., Aei=200., Aie=125., Aii=125., blocking_θI=29., stim_strength=2.), # monotonic non-prop. (not sus.)
        (Aee=120., Aei=200., Aie=125., Aii=125., blocking_θI=29., stim_strength=3.), # monotonic mid prop. (not sus.)
        (Aee=120., Aei=200., Aie=125., Aii=125., blocking_θI=29., stim_strength=20.), # monotonic prop. (not sus.)
        (Aee=120., Aei=80., Aie=125., Aii=125., blocking_θI=29., stim_strength=2.), # monotonic non-prop. (sus.)
        (Aee=120., Aei=80., Aie=125., Aii=125., blocking_θI=29., stim_strength=3.), # monotonic mid prop. (sus.)
        (Aee=120., Aei=80., Aie=125., Aii=125., blocking_θI=29., stim_strength=20.), # monotonic prop. (sus.)
        (Aee=40., Aei=240., Aie=125., Aii=125., blocking_θI=9., stim_strength=2.), # dep-block non-prop. (not sus.)
        (Aee=40., Aei=240., Aie=125., Aii=125., blocking_θI=9., stim_strength=3.), # dep-block mid prop. (not sus.)
        (Aee=40., Aei=240., Aie=125., Aii=125., blocking_θI=9., stim_strength=20.), # dep-block prop. (not sus.)
        (Aee=40., Aei=240., Aie=125., Aii=125., blocking_θI=29., stim_strength=2.), # monotonic non-prop. (not sus.)
        (Aee=40., Aei=240., Aie=125., Aii=125., blocking_θI=29., stim_strength=3.), # monotonic mid prop. (not sus.)
        (Aee=40., Aei=240., Aie=125., Aii=125., blocking_θI=29., stim_strength=20.), # monotonic prop. (not sus.)
        (Aee=300., Aei=60., Aie=125., Aii=125., blocking_θI=9., stim_strength=2.), # dep-block non-prop. (always)
        (Aee=300., Aei=60., Aie=125., Aii=125., blocking_θI=9., stim_strength=3.), # dep-block mid prop. (always)
        (Aee=300., Aei=60., Aie=125., Aii=125., blocking_θI=9., stim_strength=20.), # dep-block prop. (always)
        (Aee=300., Aei=60., Aie=125., Aii=125., blocking_θI=29., stim_strength=2.), # monotonic non-prop. (always)
        (Aee=300., Aei=60., Aie=125., Aii=125., blocking_θI=29., stim_strength=3.), # monotonic mid prop. (always)
        (Aee=300., Aei=60., Aie=125., Aii=125., blocking_θI=29., stim_strength=20.) # monotonic prop. (always)
    ];
    session_name=session_name,
    session_id=session_id
)

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