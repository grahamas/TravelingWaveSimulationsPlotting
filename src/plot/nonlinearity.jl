using Simulation73, NeuralModels
using NeuralModels: AbstractSigmoidNonlinearityParameter, AbstractDifferenceOfSigmoidsParameter

TOL = 0.001

function _bounds_fns(nonl::NeuralModels.NonlinearityParameterWrapperAction, args...)
    _bounds_fns(nonl.p, args...)
end

function _bounds_fns(nonl::AbstractSigmoidNonlinearityParameter, tol=TOL)
    left = (xs) -> findfirst(xs .> tol)
    right = (xs) -> findfirst(xs .> (1.0 - tol))
    return (left, right)
end

function _bounds_fns(nonl::Union{AbstractDifferenceOfSigmoidsParameter,AbstractDifferenceOfSigmoidsAction}, tol=TOL)
    left = (xs) -> findfirst(xs .> tol)
    right = (xs) -> findfirst(xs[left(xs):end] .< tol)
    return (left, right)
end

function _auto_range_idx!(output, nonl::AbstractNonlinearityAction)
    (left_fn, right_fn) = _bounds_fns(nonl)
    left_idx, right_idx = (left_fn(output), right_fn(output))
    @show right_idx
    # TODO error messages
    # left_idx = left_idx isa Int ? max(left_idx-20,firstindex(output)) : firstindex(output)
    left_idx = firstindex(output)
    right_idx = right_idx isa Int ? min(right_idx+((right_idx - left_idx)÷4), lastindex(output)) : lastindex(output)
    return (left_idx:right_idx)
end

function _auto_range_idx!(outputs::AbstractArray{<:AbstractArray}, ps::AbstractArray)
    all_idxs = _auto_range_idx!.(outputs, ps)
    firsts = [idx[begin] for idx in all_idxs]
    lasts = [idx[end] for idx in all_idxs]
    @show lasts
    return (minimum(firsts):maximum(lasts))
end

_potential_range(::AbstractNonlinearityParameter) = 0.:0.01:100.
function _potential_range(sig::AbstractSigmoidNonlinearityParameter)
    @assert sig.a > 0
    nearly_one = NeuralModels.inverse_simple_sigmoid_fn(0.999, sig.a, sig.θ)
    radius = nearly_one - sig.θ
    return range(sig.θ-radius, stop=nearly_one, length=100)
end
function _potential_range(dos::AbstractDifferenceOfSigmoidsParameter)
    @assert get_firing_sigmoid(dos).θ < get_blocking_sigmoid(dos).θ
    firing_range = _potential_range(get_firing_sigmoid(dos))
    blocking_range = _potential_range(get_blocking_sigmoid(dos))
    return range(firing_range[begin], stop=blocking_range[end], length=1000)
end
function _potential_range(arr::AbstractArray{<:AbstractNonlinearityParameter})
    ranges = _potential_range.(arr)
    left = minimum(first.(ranges))
    right = maximum(last.(ranges))
    return range(left, stop=right, length=1000)
end

function _auto_range(nonl_param::AbstractNonlinearityParameter, potential_range=_potential_range(nonl_param))
    output = copy(potential_range)
    dummy_space = CompactLattice(0.0, 1.0, length(potential_range))
    nonl = nonl_param(dummy_space)
    nonl(output, nothing, nothing)
    idxs = _auto_range_idx!(output, nonl)
    (potential_range[idxs], output[idxs])
end
    
function _auto_range(nonl_params::AbstractArray{<:AbstractNonlinearityParameter}, potential_range=_potential_range(nonl_params))
    @show potential_range
    @show nonl_params
    dummy_space = CompactLattice(first(potential_range), last(potential_range), length(potential_range))
    nonls = map(np -> np(dummy_space), nonl_params)
    max_outputs = [nonl(collect(potential_range), nothing, nothing) for nonl in nonls]
    idxs = _auto_range_idx!(max_outputs, nonls)
    @show idxs
    ([potential_range[idxs] for _ in max_outputs], [output[idxs] for output in max_outputs])
end    

function convert_arguments(nonl::AbstractNonlinearityParameter)
    _auto_range(nonl)
end

function AbstractPlotting.convert_arguments(P::Type{<:AbstractPlot}, pops::Simulation73.AbstractPOneD{NPOPS,<:AbstractNonlinearityParameter}) where NPOPS
    return _auto_range(Simulation73.array(pops))
end

plot_nonlinearity!(fig::Figure, sim::Simulation, args...; kwargs...) = plot_nonlinearity!(fig, sim.model, args...; kwargs...)

plot_nonlinearity!(fig::Figure, model::AbstractModel, args...; kwargs...) = plot_nonlinearity!(fig, model.nonlinearity, args...; pop_names=model.pop_names, kwargs...)

using Colors
function plot_nonlinearity!(fig::Figure, nonlinearities::Simulation73.AbstractPopulationP;
        pop_names,
        title="Population activation functions")

    nonl_xs, nonl_ys = _auto_range(nonlinearities |> Simulation73.array)
    ax = MakieLayout.Axis(fig, xlabel="input (a.u.)", ylabel="pop. activity (proportion)")
    if title !== nothing
        ax.title = title
    end
    colors = distinguishable_colors(length(pop_names), parse(Colorant, ax.attributes[:backgroundcolor][]), dropseed=true)
    plots = [lines!(ax, xs, ys, width=3, color=color) 
        for (xs, ys, color) in zip(nonl_xs, nonl_ys, colors)]
    ylims!(ax, 0.0, 1.0)
    hidespines!(ax, :l)
    ax.yticks = [0, 1]
    legend_names = ["$(pop)" for pop in pop_names]
    axislegend(ax, plots, legend_names,
                  tellheight=false, tellwidth=false,
                  position=:lt, orientation=:vertical)

    ax
end

