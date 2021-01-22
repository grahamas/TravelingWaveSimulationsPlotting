function figure_example!(scene::Scene, mdb_path::String,
        fixing_mods::NamedTuple{NAMES}) where NAMES
    prototype_name, all_mods = read_params_from_data_path(mdb_path)
    figure_example!(scene, prototype_name, all_mods, fixing_mods)
end

function figure_example!(scene::Scene, prototype_name, all_mods, fixing_mods::NamedTuple{NAMES}; title=string(fixing_mods), kwargs...) where NAMES
    fixed_mods = Dict(key => val for (key,val) in pairs(all_mods) if length(val)== 1)
    unfixed_mod_names = [key for (key, val) in pairs(all_mods) if length(val) > 1]
    @assert Set(NAMES) == Set(unfixed_mod_names) "$(Set(NAMES) - Set(unfixed_mod_names)): all varying parameters must be fixed by example, and no extraneous mods are allowed"

    # other_opts = Dict() makes sure it saves all frames
    mods = (fixed_mods..., fixing_mods..., save_idxs=nothing, save_on=true)
    #write_modifications!(plots_subdir, mods, unique_id)
    prototype = get_prototype(prototype_name)
    (mod_name, exec) = execute_single_modification(prototype, mods)

    layout = exec_heatmap!(scene, exec; colorrange=(0.0,0.5), title=title, kwargs...)
    # layout[end+1,2] = LText(scene, string(exec.simulation.global_reduction(exec.solution).propagation), tellwidth=false)

    return layout
end

function figure_example(prototype_name, all_mods, fixing_mods;
        scene_resolution=(800,400), kwargs...)
    scene, layout = layoutscene(resolution=scene_resolution)
    layout[1,1:2] = figure_example!(scene, prototype_name, all_mods, fixing_mods; kwargs...)
    fixing_mods_strs = ["$(name)=$(val)" for (name, val) in pairs(fixing_mods)]
    return (scene, layout, join(fixing_mods_strs, "_"))
end

function figure_example(mdb_path, fixing_mods; kwargs...)
    prototype_name, all_mods = read_params_from_data_path(mdb_path)
    figure_example(prototype_name, all_mods, fixing_mods; kwargs...)
end

function multiple_figure_examples(mdb_path, fixing_mods_list; kwargs...)
    prototype_name, all_mods = read_params_from_data_path(mdb_path)
    map(fixing_mods_list) do fixing_mods
        figure_example(prototype_name, all_mods, fixing_mods; kwargs...)
    end
end

function multiple_figure_examples!(scene, mdb_path, fixing_mods_list; kwargs...)
    prototype_name, all_mods = read_params_from_data_path(mdb_path)
    map(fixing_mods_list) do fixing_mods
        fixing_mods_strs = ["$(name)=$(val)" for (name, val) in pairs(fixing_mods)]
        (scene, figure_example!(scene, prototype_name, all_mods, fixing_mods; kwargs...), join(fixing_mods_strs, "_"))
    end
end

export multiple_figure_examples, multiple_figure_examples!