using TravelingWaveSimulationsPlotting, TravelingWaveSimulations
using Dates
using Makie, AbstractPlotting

session_name = "wider_strength_depthreshold_A"
session_id = "$(Dates.now())"

metasweep_A_fpath = TravelingWaveSimulations.get_recent_simulation_data_path(joinpath(homedir(), "data", "ring_normed_blocking", session_name))

clean_theme = Theme(
    fontsize=11,
    Axis = (
        #fontsize=8,
        # bottomspinevisible = false,
        # leftspinevisible = false,
        # topspinevisible = false,
        # rightspinevisible = false,
        xticksvisible = false,
        yticksvisible = false
    )
)
set_theme!(clean_theme)

kwargs = (no_labels=true,
title=nothing,
axis_kwargs = (
    xticksvisible=false,
))

function figure_propagation_susceptability_by_depblock_and_stim!(scene, fpath; kwargs...)
    layout = GridLayout()

    no_prop_param_list = [
        (Aee=120., Aei=200., Aie=125., Aii=125., blocking_θI=29., stim_strength=2.), # monotonic non-prop. (not sus.)
        (Aee=120., Aei=200., Aie=125., Aii=125., blocking_θI=29., stim_strength=3.), # monotonic mid prop. (not sus.)
        (Aee=120., Aei=200., Aie=125., Aii=125., blocking_θI=29., stim_strength=20.), # monotonic prop. (not sus.)
        (Aee=40., Aei=240., Aie=125., Aii=125., blocking_θI=9., stim_strength=2.), # dep-block non-prop. (not sus.)
        (Aee=40., Aei=240., Aie=125., Aii=125., blocking_θI=9., stim_strength=3.), # dep-block mid prop. (not sus.)
        (Aee=40., Aei=240., Aie=125., Aii=125., blocking_θI=9., stim_strength=20.) # dep-block prop. (not sus.)
    ]
    sus_prop_param_list = [
        (Aee=120., Aei=200., Aie=125., Aii=125., blocking_θI=9., stim_strength=2.), # dep-block non-prop. (sus.)
        (Aee=120., Aei=200., Aie=125., Aii=125., blocking_θI=9., stim_strength=3.), # dep-block mid prop. (sus.)
        (Aee=120., Aei=200., Aie=125., Aii=125., blocking_θI=9., stim_strength=20.), # dep-block prop. (sus.)
        (Aee=120., Aei=80., Aie=125., Aii=125., blocking_θI=29., stim_strength=2.), # monotonic non-prop. (sus.)
        (Aee=120., Aei=80., Aie=125., Aii=125., blocking_θI=29., stim_strength=3.), # monotonic mid prop. (sus.)
        (Aee=120., Aei=80., Aie=125., Aii=125., blocking_θI=29., stim_strength=20.) # monotonic prop. (sus.)
    ]

    all_prop_param_list = [
        (Aee=300., Aei=60., Aie=125., Aii=125., blocking_θI=9., stim_strength=2.), # dep-block non-prop. (always)
        (Aee=300., Aei=60., Aie=125., Aii=125., blocking_θI=9., stim_strength=3.), # dep-block mid prop. (always)
        (Aee=300., Aei=60., Aie=125., Aii=125., blocking_θI=9., stim_strength=20.), # dep-block prop. (always)
        (Aee=300., Aei=60., Aie=125., Aii=125., blocking_θI=29., stim_strength=2.), # monotonic non-prop. (always)
        (Aee=300., Aei=60., Aie=125., Aii=125., blocking_θI=29., stim_strength=3.), # monotonic mid prop. (always)
        (Aee=300., Aei=60., Aie=125., Aii=125., blocking_θI=29., stim_strength=20.) # monotonic prop. (always)
    ]

    blocking_θI_values = [9.0, 29.0]

    # stim_strength_values = [2.0, 3.0, 20.0]
    stim_strength_values = [2.0, 20.0]
    n_stim_strength_values = length(stim_strength_values)

    params_scene_layout_name_list_by_prop = Dict()
    params_scene_layout_name_list_by_prop["Not sus"] = zip(no_prop_param_list, multiple_figure_examples!(scene, fpath, no_prop_param_list; kwargs...))
    params_scene_layout_name_list_by_prop["Sus"] = zip(sus_prop_param_list, multiple_figure_examples!(scene, fpath, sus_prop_param_list; kwargs...))
    params_scene_layout_name_list_by_prop["Always"] = zip(all_prop_param_list, multiple_figure_examples!(scene, fpath, all_prop_param_list; kwargs...))
    n_prop_conditions = length(params_scene_layout_name_list_by_prop)

    for (θI_idx, θI_value) in enumerate(blocking_θI_values)
        sublayout = GridLayout()
        for (prop_cond_idx, (prop_name, params_scene_layout_name_list)) in enumerate(pairs(params_scene_layout_name_list_by_prop))
            subsublayout = map(stim_strength_values) do stim_strength
                only(_layout for 
                    (_param, (_scene, _layout, _name)) in 
                        params_scene_layout_name_list 
                    if _param[:blocking_θI] == θI_value &&
                        _param[:stim_strength] == stim_strength
                )
            end
            sublayout[prop_cond_idx,2] = LText(scene, prop_name, tellheight=false, rotation=pi/2)
            sublayout[prop_cond_idx,3:length(subsublayout)+2] = subsublayout
        end
        sublayout[1:size(sublayout)[1],1] = LText(scene, "θI = $θI_value", tellheight=false, rotation=pi/2)

        sublayout[1,2:size(sublayout)[2]-1] = [LText(scene, "S = $stim_strength", tellwidth=false) for stim_strength in stim_strength_values]
        #layout[2+(θI_idx-1)*size(sublayout)[1]:1+θI_idx*size(sublayout)[1],
        #    2:size(sublayout)[2]+1] = sublayout
        layout[2:size(sublayout)[2]+1] = sublayout
    end
    trim!(layout)

    return layout
end

plot_and_save(figure_propagation_susceptability_by_depblock_and_stim!, metasweep_A_fpath; scene_resolution=(1000, 400), kwargs...)