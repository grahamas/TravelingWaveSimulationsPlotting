function metasweep_plot!(scene, (x_axis,)::NTuple{1,<:AbstractVector{T}}, vals::AbstractVector{T};  lims = nothing,
        title="fixme", plotted_var_name::AbstractString, 
        metasweep_var_names::NTuple{1}, kwargs...) where {T<:Number}
    ax = LAxis(scene,  xlabel=metasweep_var_names[1],
        ylabel=plotted_var_name)
        vals_copy = NaN_if_missing.(Val(T), vals)
    plot!(ax, collect(x_axis), vals_copy; kwargs...)
    if lims !== nothing
        ylims!(ax, lims...)
    end
    return ax
end

function metasweep_plot!(scene, (x_axis,)::NTuple{1,<:AbstractVector{T}}, stats_and_estimates::Vector{<:Tuple{<:Any,<:Estimated}}; 
        title="", plotted_var_name::AbstractString, lims = nothing, 
        metasweep_var_names::NTuple{1}, kwargs...) where T
    # FIXME add title back
    ax = LAxis(scene, xlabel=metasweep_var_names[1], 
        ylabel=plotted_var_name)
    full_statistics = map(est -> NaN_if_missing(Val(T), est[1]), stats_and_estimates)
    est_means = map(est -> est[2].mean, stats_and_estimates)
    est_bands = map(est -> est[2].band, stats_and_estimates)
    all(isnan.(full_statistics)) && @warn "All NaN for full statistics!"
    plot!(ax, x_axis, full_statistics; kwargs...)
    errorbars!(ax, x_axis, est_means, est_bands)
    plot!(ax, x_axis, est_means, color=:red, markersize=4)
    if lims !== nothing
        ylims!(ax, lims...)
    end
    return ax
end

function metasweep_plot!(scene::Scene, (x_axis, y_axis)::Tuple{<:AbstractVector{T},<:AbstractVector{T}}, vals::AbstractMatrix{T};  lims=nothing,
        title="fixme", plotted_var_name::AbstractString, 
        metasweep_var_names::NTuple{2}, kwargs...) where {T<:Number}
    vals = collect(vals)
    layout = GridLayout()
    layout[1,1] = ax = LAxis(scene,  
        xlabel=metasweep_var_names[1],
        ylabel=metasweep_var_names[2])
    vals_copy = NaN_if_missing.(Val(T), vals)
    plt = heatmap!(ax, collect(x_axis), collect(y_axis), vals_copy; kwargs...)
    if lims !== nothing
        plt.colorrange = lims
    end
    tightlimits!(ax)
    layout[1,2] = LColorbar(scene, plt, width=20, 
        label=plotted_var_name)
    return layout
end

function metasweep_plot!(scene::Scene, (x_axis, y_axis)::Tuple{<:AbstractVector{T},<:AbstractVector{T}}, stats_and_estimates::AbstractMatrix{<:Tuple{<:Any,<:Estimated}};  lims=nothing,
        title="fixme", plotted_var_name::AbstractString, 
        metasweep_var_names::NTuple{2}, kwargs...) where {T<:Number}
    layout = GridLayout()
    layout[1,1] = ax = LAxis(scene,  
        xlabel=metasweep_var_names[1],
        ylabel=metasweep_var_names[2])
    full_statistics = map(est -> NaN_if_missing(Val(T), est[1]), stats_and_estimates)
    plt = heatmap!(ax, collect(x_axis), collect(y_axis), full_statistics; kwargs...)
    tightlimits!(ax)
    if lims !== nothing
        plt.colorrange = lims
    end
    layout[1,2] = LColorbar(scene, plt, width=20, 
        label=plotted_var_name)
    return layout
end


function extract_and_plot_metasweep!(scene, metasweep_values, bootstraps, 
        extract_fn::Function; band_fn::Function, plot_kwargs...)
    estimates = map(bootstraps) do bootstrap
        full_estimate = extract_fn(bootstrap.result_from_full)
        non_missing = skipmissing([extract_fn(result) for result in bootstrap.results_from_subsampling]) |> collect
        if length(non_missing) < 2
            return (full_estimate, Estimated(missing, missing))
        end
        estimate_band = band_fn(non_missing) 
        return (full_estimate, Estimated(mean(non_missing), estimate_band))
    end
    return metasweep_plot!(scene, metasweep_values, estimates; plot_kwargs...)
end

function plot_averaged_to_axes!(scene, data, axis_syms::NTuple{N,Symbol}) where N
    averaged_to_axes = _collapse_to_axes(data, axis_syms...)
    axis_values = map(sym -> get_axis_values(data, sym), axis_syms)
    layout = metasweep_plot!(scene, axis_values, averaged_to_axes; plotted_var_name="average", metasweep_var_names=string.(axis_syms))
    return layout
end

function plot_averaged_to_axes(data, axes::Vector{NTuple{N,Symbol}}; scene_resolution=(800*length(axes),600)) where N
    scene, layout = layoutscene(resolution=scene_resolution)
    layout[1,1:N] = map(axis_syms -> plot_averaged_to_axes!(scene, data, axis_syms), axes)
    return (scene, layout)
end

plot_averaged_to_axes(data, axes::NTuple) = plot_averaged_to_axes(data, [axes])

function multiple_averaged_to_axes(data::AbstractArray; axes, multi_param_syms, kwargs...)
    multi_param_values = map(sym -> get_axis_values(data, sym), multi_param_syms)
    multi_param_coords = product(multi_param_values...)
    scene_layout_name_list = map(multi_param_coords) do coord
        named_coord = NamedTuple(name => val for (name, val) in zip(multi_param_syms, coord))
        data_slice = getindex(data; named_coord...)
        scene, layout = plot_averaged_to_axes(data_slice, axes; kwargs...)
        coord_strings = ["$(sym)=$(val)" for (sym,val) in pairs(named_coord)]
        return (scene, layout, join(coord_strings, "_"))
    end
    return scene_layout_name_list
end

multiple_averaged_to_axes(mdb_path::String; property_sym, kwargs...) = multiple_averaged_to_axes(TravelingWaveSimulations.load_classifications(mdb_path)[property_sym]; kwargs...)

export multiple_averaged_to_axes
