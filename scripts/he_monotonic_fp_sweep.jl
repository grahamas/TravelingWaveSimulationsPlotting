using TravelingWaveSimulationsPlotting, TravelingWaveSimulations
using WilsonCowanModel, Simulation73
using Dates
using Contour
using AxisIndices, IterTools
using Makie

A_range = 0.1:0.1:1.5
static_mods = (
    α=(0.4, 0.7),
    n_lattice = 2,
    save_idxs=nothing, save_on=true, saveat=0.1
) 
sweep_calculate_fixedpoints_and_plot(; 
    nonl_type="monotonic", 
    A_range=A_range, 
    static_mods=static_mods
)