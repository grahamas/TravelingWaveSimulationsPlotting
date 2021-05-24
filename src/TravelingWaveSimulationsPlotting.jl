module TravelingWaveSimulationsPlotting

using DrWatson

using Makie
using Makie.MakieLayout
using TravelingWaveSimulations, Simulation73, NeuralModels, Simulation73Plotting, WilsonCowanModel
using Interpolations, DiffEqOperators, Optim, LinearAlgebra
using AxisIndices
using StaticArrays
using Contour: contour, lines, coordinates, Curve2
using IterTools: product
using LoopVectorization
using LinearAlgebra
using ForwardDiff
using OffsetArrays
using Roots
using Dates
using Colors

using DataStructures: MutableLinkedList, ListNode, length

using NeuralModels: AbstractSigmoidNonlinearityParameter, AbstractDifferenceOfSigmoidsParameter

using WilsonCowanModel: wcm_du_defn, wcm_dv_defn, WCMParams, AbstractNullclineParams, get_nullcline_params
using Makie: @lift

include("util/axisarray.jl")
include("util/lines.jl")
#include("util/bootstrap.jl")
include("util/fitting_sigmoid.jl")
export PointVectorLine, Bootstrapped, Estimated
using Statistics: std, mean
using HypothesisTests: OneSampleTTest, confint
include("util/sweep.jl")
export wcm_sweep_calculate_fixedpoints
include("util/metrics.jl")
export epilepsy_metric
include("util/fixedpoints.jl")
export calculate_fixedpoints,
    calculate_fixedpoints!
include("util/fixedpoint_stability.jl")

include("space_reduction.jl")
include("phase_space.jl")

include("plotting.jl")
export heatmap_sweep_with_target,
    figure_contrast_monotonic_blocking_all, 
    figure_example_contrast_monotonic_blocking_all,
    save_figure_example_contrast_monotonic_blocking_all,
    plot_and_save, layout_plot, figure_plot
include("plot/axisarray.jl")
export axisarray_heatmap!
include("plot/sweep.jl")
export sweep_calculate_fixedpoints_and_plot
include("plot/nullclines.jl")
export plot_nullclines!
include("plot/connectivity.jl")
export plot_connectivity!
include("plot/nonlinearity.jl")
export plot_nonlinearity!
include("plot/stimulus.jl")
export plot_stimulus!

include("figures/binarization.jl")
include("figures/sigmoid_fit.jl")
export figure_metasweep_sigmoid_fits!,
    figure_metasweep_binarization_counts!
include("figures/example.jl")
export figure_example!

end # module
