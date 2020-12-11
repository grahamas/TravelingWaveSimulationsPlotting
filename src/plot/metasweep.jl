function metasweep_plot!(scene, (x_axis,)::NTuple{1,<:AbstractVector{T}}, vals::AbstractVector{T}; 
        title="fixme", plotted_var_name::AbstractString, 
        metasweep_var_names::NTuple{1}, kwargs...) where {T<:Number}
    ax = LAxis(scene,  xlabel=metasweep_var_names[1],
        ylabel=plotted_var_name)
        vals_copy = NaN_if_missing.(Val(T), vals)
    plot!(ax, collect(x_axis), vals_copy; kwargs...)
    return ax
end

function metasweep_plot!(scene, (x_axis,)::NTuple{1,<:AbstractVector{T}}, stats_and_estimates::Vector{<:Tuple{<:Any,<:Estimated}}; 
        title="", plotted_var_name::AbstractString, 
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
    return ax
end

function metasweep_plot!(scene::Scene, (x_axis, y_axis)::Tuple{<:AbstractVector{T},<:AbstractVector{T}}, vals::Matrix{T};  
        title="fixme", plotted_var_name::AbstractString, 
        metasweep_var_names::NTuple{2}, kwargs...) where {T<:Number}
    layout = GridLayout()
    layout[1,1] = ax = LAxis(scene,  
        xlabel=metasweep_var_names[1],
        ylabel=metasweep_var_names[2])
    vals_copy = NaN_if_missing.(Val(T), vals)
    plt = heatmap!(ax, collect(x_axis), collect(y_axis), vals_copy; kwargs...)
    tightlimits!(ax)
    layout[1,2] = LColorbar(scene, plt, width=20, 
        label=plotted_var_name)
    return layout
end

function metasweep_plot!(scene::Scene, (x_axis, y_axis)::Tuple{<:AbstractVector{T},<:AbstractVector{T}}, stats_and_estimates::AbstractMatrix{<:Tuple{<:Any,<:Estimated}};  
        title="fixme", plotted_var_name::AbstractString, 
        metasweep_var_names::NTuple{2}, kwargs...) where {T<:Number}
    layout = GridLayout()
    layout[1,1] = ax = LAxis(scene,  
        xlabel=metasweep_var_names[1],
        ylabel=metasweep_var_names[2])
    full_statistics = map(est -> NaN_if_missing(Val(T), est[1]), stats_and_estimates)
    plt = heatmap!(ax, collect(x_axis), collect(y_axis), full_statistics; kwargs...)
    tightlimits!(ax)
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