using Dates

function save_metasweep(mdb_path, args...; 
        unique_id="$(Dates.now())",
        plots_subdir=plotsdir("metasweep_$(unique_id)"),
        kwargs...)
    scene, _ = figure_metasweep(mdb_path, args...; 
                                      plots_subdir=plots_subdir,
                                      kwargs...)
    fname = "metasweep_$(unique_id).png"
    mkpath(plotsdir(plots_subdir))
    save_name = plotsdir(plots_subdir, fname)
    @info "saving $save_name"
    Makie.save(save_name, scene)
    return scene
end

function figure_metasweep(args...; 
                        scene_resolution=(600, 600), kwargs...)
    scene, layout = layoutscene(resolution=scene_resolution)

    layout[1,1] = plot_metasweep!(scene, args...; kwargs...)

    return (scene, layout)
end

function metasweep_plot!(scene, (x_axis,)::NTuple{1}, args...; title="", kwargs...)
    ax = LAxis(scene, title=title)
    plot!(ax, x_axis, args...; kwargs...)
    return ax
end

function metasweep_plot!(scene::Scene, (x_axis, y_axis)::NTuple{2}, args...; title="", kwargs...)
    layout = GridLayout()
    layout[1,1] = ax = LAxis(scene, title=title)
    plt = heatmap!(ax, x_axis, y_axis, args...; kwargs...)
    tightlimits!(ax)
    layout[1,2] = LColorbar(scene, plt, width=20)
    return layout
end

function plot_metasweep!(scene::Scene, mdb_path::String, 
        (x_sym, y_sym)::Tuple{Symbol,Symbol}, 
        metasweep_syms::NTuple{N,Symbol}, 
        property_sym::Symbol; plots_subdir=nothing) where N
    data = TravelingWaveSimulations.load_ExecutionClassifications(AbstractArray, mdb_path)[property_sym]

    metasweep_values = map(sym -> axes_keys(data)[dim(data, sym)], metasweep_syms)
    metasweep_value_pairs = product(metasweep_values...)
    getindex_metasweep_dim(arr, idxs) = getindex(arr; Dict(sym => idx for (sym, idx) in zip(metasweep_syms, idxs))...)

    (xs, ys) = axes_keys(_collapse_to_axes(getindex_metasweep_dim(data, first(metasweep_value_pairs)), x_sym, y_sym))

    fitted_sigmoids_and_threshold_locs = map(metasweep_value_pairs) do metasweep_value_pair
        metasweep_slice = getindex_metasweep_dim(data, metasweep_value_pair)
        # flattened_metasweep_slice = _collapse_to_axes(metasweep_slice, x_sym, y_sym)
        # if plots_subdir !== nothing
            # save_reduce_2d_and_steepest_line_and_histogram((x_sym, y_sym), flattened_metasweep_slice,
                    # property_sym, "$metasweep_value_pair", plots_subdir; facet_title="$metasweep_value_pair")
        # end
        phase_space_fitted_sigmoid = bootstrap_sigmoid_fit(metasweep_slice, x_sym, y_sym) do subsampled0_slice
            flattened_subsampled_slice = _collapse_to_axes(subsampled_slice, x_sym, y_sym)
            dists, vals, locs, lin = reduce_normal_to_halfmax_contour(flattened_subsampled_slice, slice)
            fitted_sigmoid = fit_sigmoid(vals, dists)
            phase_space_sigmoid = SigmoidOnSlice(fitted_sigmoid, lin)
        end
        return phase_space_fitted_sigmoid
    end
    fitted_sigmoids = [a[1] for a in fitted_sigmoids_and_threshold_locs]
    threshold_locs = [a[2] for a in fitted_sigmoids_and_threshold_locs]

    ifnonmissing(x)::Float64 = ismissing(x) ? NaN64 : x
    changes = map(x -> x.change |> ifnonmissing, fitted_sigmoids)
    if all(isnan.(changes))
        error("no sigmoids fit")
    end
    slopes = map(x -> x.slope |> ifnonmissing, fitted_sigmoids)
    thresholds = map(x -> x.threshold |> ifnonmissing, fitted_sigmoids)
    x_threshold_locs = [ismissing(a) ? NaN : a[1] for a in threshold_locs]
    y_threshold_locs = [ismissing(a) ? NaN : a[2] for a in threshold_locs]
    errors = map(x -> x.error |> ifnonmissing, fitted_sigmoids)

    set_theme!(LAxis=(textsize=5,), LText=(tellwidth=false, tellheight=false))
    layout = GridLayout(resolution=(1200,1200))
    layout[1,1] = changes_layout = metasweep_plot!(scene, metasweep_values, changes, title="boundaries delta")
    layout[1,2] = slopes_layout = metasweep_plot!(scene, metasweep_values, slopes, title="slopes")
    layout[2,2] = errors_layout = metasweep_plot!(scene, metasweep_values, errors, title="errors")
    layout[3,1] = x_thresholds_layout = metasweep_plot!(scene, metasweep_values, x_threshold_locs, title="$x_sym θ")
    layout[3,2] = y_thresholds_layout = metasweep_plot!(scene, metasweep_values, y_threshold_locs, title="$y_sym θ")

    return layout
end

# offset - A sigmoid(-a (x - theta))
# where x is distance along PointVectorLine(θ, φ)
struct SigmoidOnSlice{T}
    offset::T
    A::T
    a::T
    θ::NTuple{T,2}
    φ::NTuple{T,2}
    error::T
end
function SigmoidOnSlice(sig::FittedSigmoid{T}, lin::PointVectorLine)
    SigmoidOnSlice{T}(
        sig.left_val,
        sig.change,
        sig.slope,
        point_from_distance(lin, sig.threshold),
        Tuple(lin.vector),
        sig.error
    )
end










function heatmap_sweep_with_target(sweep::AbstractArray,
			target_mods_nt::NamedTuple{mod_names}, 
			prototype_name,
            sim_name;
            fixed_mods=Dict(),
            title,
			plot_color=:magma,
    		plot_side_size = 350 * (length(mod_names) - 1)
		) where {mod_names}
    mod_values = keys.(axes(sweep))
    mod_names_str = [string(name) for name in mod_names]
    line_prototype = get_prototype(prototype_name)
    all_dims = 1:length(mod_names)
    slices_2d = IterTools.subsets(all_dims, Val{2}())
    plot_size = (plot_side_size, plot_side_size)
    scene, layout = layoutscene(resolution=plot_size)

    heatmaps = map(slices_2d) do (x,y)
        (x,y) = x < y ? (x,y) : (y,x)
        my = mod_values[y] |> collect
        mx = mod_values[x] |> collect
        collapsed_dims = Tuple(setdiff(all_dims, (x,y)))
        sweep_2d_mean = Simulation73.avg_across_dims(sweep, collapsed_dims)
        
        @assert size(sweep_2d_mean) == length.((mod_values[x], mod_values[y]))
        
        layout[x,y] = ax = LAxis(scene); 
        tightlimits!(ax)
        
        heatmap = Makie.heatmap!(ax, my, mx, sweep_2d_mean', colorrange=(0,1))
        
        Makie.scatter!(ax, [target_mods_nt[y]], [target_mods_nt[x]], color=:red, markersize=5)
        
        heatmap
    end
    savedir(bn) = plotsdir(prototype_name, sim_name, mods_filename(; fixed_mods..., target_mods_nt...), bn)
    mkpath(savedir("") |> dirname)
    
    layout[:,1] = LText.(scene, mod_names_str[1:end-1], tellheight=false, rotation=pi/2)
    layout[end+1,2:end] = LText.(scene, mod_names_str[2:end], tellwidth=false)
    layout[0, :] = LText(scene, title, textsize = 30)
    cbar = layout[2:end-1, end+1] = LColorbar(scene, heatmaps[1], label = "Proportion")
    cbar.width = 25
    path = savedir("slices.png")
    Makie.save(path, scene)
    
    these_mods =  (save_idxs=nothing, other_opts=Dict(), fixed_mods..., target_mods_nt...)
    @show these_mods
    prototype = get_prototype(prototype_name)
    (mod_name, exec) = execute_single_modification(prototype, these_mods)
    #exec = execute(line_prototype(;mods..., other_opts=Dict()))
    #wp = ExecutionClassifications(exec.solution)
    #@show wp.has_propagation
    @warn "Not validating has_propagation"
    anim_path = savedir("sim_animation.mp4")
    animate_execution(anim_path, exec);
    
    scene, layout = Simulation73.exec_heatmap_slices(exec, 5, (1400,1000))
    layout[0,:] = LText(scene, "Simulation: $(target_mods_nt)\n Inhibition blocking threshold: $(these_mods[:blocking_θI])") 
    sim_heatmap_path = savedir("sim_heatmap.png")
    Makie.save( sim_heatmap_path, scene)
    
    scene, layout = Simulation73.exec_heatmap(exec)
    layout[0,:] = LText(scene, "Simulation: $(target_mods_nt)\n Inhibition blocking threshold: $(these_mods[:blocking_θI])") 
    sim_heatmap_path = savedir("sim_E_heatmap.png")
    Makie.save( sim_heatmap_path, scene)

end

function figure_contrast_monotonic_blocking_all((x_sym, y_sym)::Tuple{Symbol,Symbol}, (other_x_sym, other_y_sym)::Tuple{Symbol,Symbol},
                                     monotonic_fpath::AbstractString, 
                                     blocking_fpath::AbstractString,
                                     property_sym::Symbol;
                                     scene_resolution=(1200,1200))    
    scene, layout = layoutscene(resolution=scene_resolution)


    layout[1:3,1] = monotonic_sweep_ax = reduce_2d_and_steepest_line_and_histogram!(
                                                            scene, 
                                                            (x_sym, y_sym),
                                                            monotonic_fpath,
                                                            property_sym; 
                                                            facet_title="Monotonic")
    layout[1:3,2] = blocking_sweep_ax = reduce_2d_and_steepest_line_and_histogram!(
                                                            scene, 
                                                            (x_sym, y_sym),
                                                            blocking_fpath,
                                                            property_sym; 
                                                            facet_title="Blocking", 
                                                            hide_y=true)
    layout[1:3,3] = other_monotonic_sweep_ax = reduce_2d_and_steepest_line_and_histogram!(
                                                           scene, 
                                                           (other_x_sym, other_y_sym),
                                                           monotonic_fpath,
                                                           property_sym; 
                                                           facet_title="Monotonic")
    layout[1:3,4] = other_blocking_sweep_ax = reduce_2d_and_steepest_line_and_histogram!(
                                                          scene, 
                                                          (other_x_sym, other_y_sym),
                                                          blocking_fpath,
                                                          property_sym; 
                                                          facet_title="Blocking",
                                                          hide_y=true,
                                                          colorbar_width=25)
    layout[end+1, 1] = LText(scene, "($(x_sym), $(y_sym))", tellwidth=false)
    layout[end+1, 2] = LText(scene, "($(x_sym), $(y_sym))", tellwidth=false)
    layout[end, 3] = LText(scene, "($(other_x_sym), $(other_y_sym))", tellwidth=false)
    layout[end, 4] = LText(scene, "($(other_x_sym), $(other_y_sym))", tellwidth=false)

    layout[end-1:end-2, 0] = LText(scene, "prop. sims", rotation=pi/2, tellheight=false)

    return scene, layout
    
end

function figure_example_contrast_monotonic_blocking_all((x_sym, y_sym)::Tuple{Symbol,Symbol}, 
                                                        other_syms::Tuple{Symbol,Symbol},
                                                        example_specs::Array{<:NamedTuple},
                                     monotonic_fpath::AbstractString, 
                                     blocking_fpath::AbstractString,
                                     property_sym::Symbol; kwargs...)
    monotonic_prototype_name, monotonic_spec = read_params_from_data_path(monotonic_fpath)
    blocking_prototype_name, blocking_spec = read_params_from_data_path(blocking_fpath)
    @show blocking_spec

    @show monotonic_spec
    monotonic_prototype, blocking_prototype = get_prototype.((monotonic_prototype_name, blocking_prototype_name))

    scene_height = 450 * 3# (2 + length(example_specs))
    scene_width = 450 * 8
    scene, layout = figure_contrast_monotonic_blocking_all((x_sym, y_sym), other_syms, monotonic_fpath, blocking_fpath, property_sym; scene_resolution=(scene_width, scene_height), kwargs...)
    scene.attributes.attributes[:fontsize] = 20

    @show merge(blocking_spec, pairs(example_specs[2]))
    examples_layout = GridLayout()
    count = 0
    for spec in example_specs
        count += 1
        _, monotonic_example_exec = execute_single_modification(monotonic_prototype, merge(monotonic_spec, pairs(spec), Dict(:save_everystep => true)))
        @show ExecutionClassifications(monotonic_example_exec, n_traveling_frames_threshold=60).has_propagation
        _, blocking_example_exec = execute_single_modification(blocking_prototype, merge(blocking_spec, pairs(spec), Dict(:save_everystep => true)))
        @show ExecutionClassifications(blocking_example_exec, n_traveling_frames_threshold=60).has_propagation
        
        examples_layout[count,1] = monotonic_ax = exec_heatmap!(scene, monotonic_example_exec; clims=(0.0,0.5), no_labels=true)
        examples_layout[count,2] = blocking_ax = exec_heatmap!(scene, blocking_example_exec;
                                                           clims=(0.0,0.5), no_labels=true)
        examples_layout[count,3] = spec_text = LText(scene, 
                                                     join(["$k = $v" for (k,v) in pairs(spec)], "\n"), tellheight=false)
    end

    examples_layout[0,1] = LText(scene, "Monotonic", tellwidth=false)
    examples_layout[1,2] = LText(scene, "Blocking", tellwidth=false)
    examples_layout[end,0] = LText(scene, "space (μm)", rotation=pi/2, tellheight=false)
    examples_layout[end+1,2] = LText(scene, "time (ms)", tellwidth=false)
    
    ex_x, ex_y = size(examples_layout)
    @show size(examples_layout)
    next_col = size(layout)[2] + 1
    layout[1:4,next_col:next_col+3] = examples_layout
    
	label_a = layout[1, 2, TopLeft()] = LText(scene, "A", textsize = 35,
    	font = "Noto Sans Bold", halign = :right)
	label_b = layout[1, 4, TopLeft()] = LText(scene, "B", textsize = 35,
    	font = "Noto Sans Bold", halign = :right)
	label_c = layout[1, 6, TopLeft()] = LText(scene, "C", textsize = 35,
    	font = "Noto Sans Bold", halign = :right)

    return scene, layout
end



function save_figure_example_contrast_monotonic_blocking_all((x_sym, y_sym)::Tuple{Symbol,Symbol}, other_syms::Tuple{Symbol,Symbol}, 
                                                         example_specs::Array{<:NamedTuple},
                                     monotonic_fpath::AbstractString, 
                                     blocking_fpath::AbstractString,
                                     property_sym::Symbol,
                                     unique_id::String=""; kwargs...)
    scene, _ = figure_example_contrast_monotonic_blocking_all((x_sym, y_sym), other_syms,
                                               example_specs,
                                               monotonic_fpath,
                                               blocking_fpath,
                                               property_sym; kwargs...)
    fname = "figure_examples_contrast_monotonic_blocking_all_$(x_sym)_$(y_sym)_$(property_sym).png"
    mkpath(plotsdir(unique_id))
    @info "saving $(plotsdir(unique_id,fname))"
    Makie.save(plotsdir(unique_id,fname), scene)
end


