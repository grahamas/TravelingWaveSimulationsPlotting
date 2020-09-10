module TravelingWaveSimulationsPlotting

using DrWatson

using MakieLayout, Makie
using TravelingWaveSimulations, Simulation73, NeuralModels, Simulation73Plotting
using Interpolations, DiffEqOperators, Optim, LinearAlgebra
using AxisIndices
using StaticArrays

include("util/axisarray.jl")
include("util/lines.jl")
export PointVectorLine

include("space_reduction.jl")
include("fitting_sigmoid.jl")
include("phase_space.jl")

include("plotting.jl")
export heatmap_sweep_with_target, axisarray_heatmap!,
    figure_contrast_monotonic_blocking_all, 
    figure_example_contrast_monotonic_blocking_all,
    save_figure_example_contrast_monotonic_blocking_all,
    save_sweep_threshold
include("plot_axisarray.jl")

end # module
