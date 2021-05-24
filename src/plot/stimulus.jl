
function plot_stimulus!(fig::Figure, simulation::Simulation{T};
        ts = nothing,
        clims = nothing,
        ylabel="space (μm)",
        xlabel="time (ms)",
        title="Stimulus",
        cbar_width = 25
    ) where T
    pop_names = simulation.model.pop_names
    xs = coordinate_axes(Simulation73.reduced_space(simulation))[1] |> collect
    stimuli = Simulation73.array(simulation.model.stimulus)
    stim_actions! = [stim(simulation.space) for stim in stimuli]

    # calculate ts (times) if not provided
    ts = if ts !== nothing
        ts
    else
        stop_time = maximum(stimulus_duration.(stimuli))
        range(0, stop=stop_time, length=100)
    end

    # preallocate stimulus array
    stim_arrs = [zeros(T, length(ts), length(xs)) for stim in stimuli]

    # calculate stimulus inplace
    for (stim_arr, stim_action!) in zip(stim_arrs, stim_actions!)
        for (i_t, t) in enumerate(ts)
            @views stim_action!(stim_arr[i_t,:], nothing, t)
        end
    end
    max_val, min_val = maximum(maximum.(stim_arrs)), minimum(minimum.(stim_arrs))
    clims = if clims !== nothing
        @assert max_val <= clims[2] && min_val >= clims[1]
        clims
    else
        [min_val, max_val]
    end


    layout = GridLayout()
    axs = layout[:v] = [Makie.Axis(fig) for stim in stimuli]
    hidedecorations!.(axs[begin:end-1])
    axs[end].xlabel = xlabel
    for (i_ax, ax) in enumerate(axs)
        ax.ylabel = "$(pop_names[i_ax]): $ylabel"
    end
    axs[end].xticks = [ts[begin], ts[end]]
    axs[begin].title = title
    hs = [heatmap!(ax, ts, xs, stim_arr) for (ax, stim_arr) in zip(axs,stim_arrs)]
    tightlimits!.(axs)
    cbar = layout[length(pop_names),2] = 
            Colorbar(fig, hs[end],
                     width = cbar_width, limits = clims, 
                     ticks = clims)

    layout
end