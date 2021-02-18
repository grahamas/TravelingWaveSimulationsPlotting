using Colors

# FIXME this should really dispatch to 
#   (::Figure, ::AbstractModel, ::AbstractSpace)
#   -> (::Figure, ::SpecificConnectivity, ::AbstractSpace)
#   in case of weird connectivities
# For simplicity, I'm just assuming a kernel interface here 
function plot_connectivity!(fig::Figure, simulation::Simulation{T,M}; 
        title = "Connectivity", 
        xlabel = "distance (μm)", ylabel = "connectivity strength (a.u.)",
        xlims = nothing, ylims = nothing,
        kwargs...) where {T, M<:AbstractWilsonCowanModel{T,1}}
    ax = AbstractPlotting.Axis(fig,
              xlabel = xlabel, # FIXME encode units in simulation
              ylabel = ylabel)
    if title !== nothing
        ax.title = title
    end
    pop_names = simulation.model.pop_names
    n_pops = length(pop_names)
    colors = distinguishable_colors(n_pops,
                                    parse(Colorant, ax.attributes[:backgroundcolor][]),
                                    dropseed=true)
    raw_space = simulation.space
    midpoint_idx = NeuralModels.fft_center_idx(raw_space)
    midpoint = Simulation73.coordinates(raw_space)[midpoint_idx]
    dists = only.(differences(raw_space, midpoint))
    dists[CartesianIndex(1):midpoint_idx] .*= -1  # ATTN would be wrong if fft_center_idx were to the right

    x_idxs = if xlims !== nothing
        xlims!(ax, xlims)
        findfirst(dists .>= xlims[1]):findlast(dists .<= xlims[2])
    else
        Colon()
    end
    ylims !== nothing && ylims!(ax, ylims)
    onto_pop_idx = 1
    connectivity_array = Simulation73.array(simulation.model.connectivity)
    plots = map(pairs(IndexCartesian(), connectivity_array) |> collect) do (idx, conn)
        kern = NeuralModels.kernel(conn, raw_space)
        lines!(ax, dists[x_idxs], kern[x_idxs], color=colors[Tuple(idx)[onto_pop_idx]])
    end
    legend_names = ["onto $pop" for pop in pop_names]
    AbstractPlotting.axislegend(ax, plots[1:n_pops], legend_names,
                  #tellheight=false, tellwidth=false,
                  position=:rt, orientation=:vertical)
    tightlimits!(ax)
    return ax    
end