using TravelingWaveSimulationsPlotting, TravelingWaveSimulations
using WilsonCowanModel, Simulation73
using Dates
using Contour
using AxisIndices, IterTools
using Makie

# 7fp values:
# Aee = 1.
# Aei = 0.8
# Aie = 1.
# Aii = 0.25

A_range = 0.1:0.3:1.5
static_mods = (
    α=(0.4, 0.7), 
    firing_θI=0.2, blocking_θI=0.5, 
    n_lattice = 2,
    save_idxs=nothing, save_on=true, saveat=0.1
) 
sweep_calculate_fixedpoints_and_plot(; 
    nonl_type="blocking", 
    A_range=A_range, 
    static_mods=static_mods
)