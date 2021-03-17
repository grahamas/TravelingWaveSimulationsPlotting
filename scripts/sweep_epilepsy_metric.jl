using TravelingWaveSimulations, TravelingWaveSimulationsPlotting
using Dates

# loads blocking_fp_arr and monotonic_fp_arr
include(scriptsdir("load/he_sweep_arrs.jl"))

let fp_arrs = [("blocking", blocking_fp_arr), 
                ("monotonic", monotonic_fp_arr)],
    session_name = "he_fp_sweeps",
    session_id = "$(Dates.now())",
    colorbar_width = 5
;




end #let fp_arrs