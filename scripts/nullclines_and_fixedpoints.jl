using TravelingWaveSimulationsPlotting, TravelingWaveSimulations
using WilsonCowanModel, Simulation73
using Dates
using Contour
using AxisIndices

session_name = "nullclines_and_fixedpoints"
session_id = "$(Dates.now())"

metasweep_A_fpath = TravelingWaveSimulations.get_recent_simulation_data_path(joinpath(homedir(), "data", "ring_normed_blocking", "wider_strength_depthreshold_A"))


@time fp_arr, cmods, proto = calculate_fixedpoints(metasweep_A_fpath, (stim_strength=2.0, blocking_θI = 9.
    #, Aii=125., Aie=125.
    #, Aii=1.0, Aie=15.0
), 0.01);
