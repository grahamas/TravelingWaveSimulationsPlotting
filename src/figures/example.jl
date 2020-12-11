function figure_example!(scene::Scene, mdb_path::String,
        fixing_mods::NamedTuple{NAMES}; plots_subdir, unique_id="", animate=true) where NAMES
    prototype_name, all_mods = read_params_from_data_path(mdb_path)
    @show all_mods
    fixed_mods = Dict(key => val for (key,val) in pairs(all_mods) if length(val)== 1)
    @show fixed_mods
    unfixed_mod_names = [key for (key, val) in pairs(all_mods) if length(val) > 1]
    @assert Set(NAMES) == Set(unfixed_mod_names) "$(Set(NAMES) - Set(unfixed_mod_names)): all varying parameters must be fixed by example, and no extraneous mods are allowed"

    # other_opts = Dict() makes sure it saves all frames
    mods = (fixed_mods..., fixing_mods..., save_idxs=nothing, save_on=true)
    @show mods
    write_modifications!(plots_subdir, mods, unique_id)
    prototype = get_prototype(prototype_name)
    (mod_name, exec) = execute_single_modification(prototype, mods)
    @show mod_name

    layout = exec_heatmap!(scene, exec; clims=(0.0,0.5), title=string(fixing_mods))
    layout[end+1,2] = LText(scene, string(exec.simulation.global_reduction(exec.solution).propagation), tellwidth=false)

    if animate
        animate_execution(joinpath(plots_subdir, "$(unique_id == "" ? "" : "$(unique_id)_")animation.mp4"), exec)
    end


    return layout
end

