module TravelingWaveSimulationsPlotting

using DrWatson

using Makie
using AbstractPlotting: MakieLayout
using AbstractPlotting.MakieLayout
using TravelingWaveSimulations, Simulation73, NeuralModels, Simulation73Plotting
using Interpolations, DiffEqOperators, Optim, LinearAlgebra
using AxisIndices
using StaticArrays
using Contour: contour, lines, coordinates
using IterTools: product

include("util/axisarray.jl")
include("util/lines.jl")
#include("util/bootstrap.jl")
include("util/fitting_sigmoid.jl")
export PointVectorLine, Bootstrapped, Estimated
#using StatsBase: sample
using Statistics: std, mean
using HypothesisTests: OneSampleTTest, confint

include("space_reduction.jl")
include("phase_space.jl")

include("plotting.jl")
export heatmap_sweep_with_target,
    figure_contrast_monotonic_blocking_all, 
    figure_example_contrast_monotonic_blocking_all,
    save_figure_example_contrast_monotonic_blocking_all,
    plot_and_save, layout_plot
include("plot/axisarray.jl")
export axisarray_heatmap!
include("plot/metasweep.jl")
include("plot/nullclines.jl")
export nullclines, nullclines!,
    calculate_fixedpoints

include("figures/binarization.jl")
include("figures/sigmoid_fit.jl")
export figure_metasweep_sigmoid_fits!,
    figure_metasweep_binarization_counts!
include("figures/example.jl")
export figure_example!

end # module
