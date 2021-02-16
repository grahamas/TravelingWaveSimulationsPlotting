function metasweep_plot!(scene, (x_axis,)::NTuple{1,<:AbstractVector{T}}, vals::AbstractVector{T};  lims = nothing,
        title="fixme", plotted_var_name::AbstractString, 
        metasweep_var_names::NTuple{1}, kwargs...) where {T<:Number}
    ax = MakieLayout.Axis(scene,  xlabel=metasweep_var_names[1],
        ylabel=plotted_var_name)
        vals_copy = NaN_if_missing.(Val(T), vals)
    plot!(ax, collect(x_axis), vals_copy; kwargs...)
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
    layout[1,1] = ax = MakieLayout.Axis(scene,  
        xlabel=metasweep_var_names[1],
        ylabel=metasweep_var_names[2])
    vals_copy = NaN_if_missing.(Val(T), vals)
    plt = heatmap!(ax, collect(x_axis), collect(y_axis), vals_copy; kwargs...)
    if lims !== nothing
        plt.colorrange = lims
    end
    tightlimits!(ax)
    layout[1,2] = Colorbar(scene, plt, width=20, 
        label=plotted_var_name)
    return layout
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
