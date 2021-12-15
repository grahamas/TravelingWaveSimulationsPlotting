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

function multiple_averaged_to_axes(data::NamedDimsArray{Names}, data_dims::NamedTuple{Names}; axes, multi_param_syms, kwargs...) where Names
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

# FIXME should separate plotting part
function sweep_calculate_fixedpoints_and_plot(; 
    nonl_type,
    A_range = 0.1:0.1:1.5,
    sweeping_mods = (#@warn "Sweeping connectivity over A_range=$(A_range)";
        (Aee=A_range, Aei=A_range, Aie=A_range, Aii=A_range)
    ),
    static_mods = ( # some defaults to make fixed points easier
        n_lattice = 2,
        save_idxs=nothing, save_on=true, saveat=0.1 # this line unnecessary?
    ),
    session_name = "he_$(nonl_type)_sweep_fp",
    session_id = "$(Dates.now())")

    #@warn "Checking $(length(A_range)^4) parameterizations..."

    fp_arr = sweep_calculate_fixedpoints(
        "full_dynamics_$(nonl_type)", 
        static_mods,
        (Aee=A_range, Aei=A_range, Aie=A_range, Aii=A_range)
        ; 
        dx = 0.01
    )

    @show fp_arr[begin]
    fp_count_arr = length.(fp_arr)
    @show extrema(fp_count_arr)

    n_fps = 0:7
    fp_count = [count(fp_count_arr .== x) for x in n_fps]
    log_fp_count = log10.(fp_count)
    max_log = log_fp_count |> maximum |> mx -> ceil(Int, mx)

    sc, ly = layoutscene(); ly[1,1] = ax = MakieLayout.Axis(sc);
    barplot!(ax, n_fps, log_fp_count)
    tightlimits!(ax)
    xlims!(ax, n_fps[begin]-0.5,n_fps[end]+0.5)
    ylims!(ax, 0, max_log)
    ax.xticks=n_fps
    ax.yticks=0:max_log
    ax.ytickformat = xs -> [x > 0 ? "10^$(Int(x))" : "$x" for x in xs]
    ax.title = "Varying connectivity of $(nonl_type) WCM"
    ax.xlabel = "# fixed points"
    ax.ylabel = "# models"
    hidespines!(ax, :t, :r)
    hidedecorations!(ax, label=false, ticklabels=false, ticks=false)
    mkpath(plotsdir("$(session_name)_$(session_id)"))
    save(plotsdir("$(session_name)_$(session_id)", "he_$(nonl_type)_fp.png"), sc)
    
    fp_arr
end