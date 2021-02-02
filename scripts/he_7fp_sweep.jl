using TravelingWaveSimulationsPlotting, TravelingWaveSimulations
using WilsonCowanModel, Simulation73
using Dates
using Contour
using AxisIndices, IterTools

@show "running"

fp_arr = let session_name = "he_7fp_sweep",
    session_id = "$(Dates.now())";

# 7fp values:
# Aee = 1.
# Aei = 0.8
# Aie = 1.
# Aii = 0.25
A_range = 0.1:0.1:1.5
@warn "Checking $(length(A_range)^4) parameterizations..."

prototype = get_prototype("full_dynamics_blocking")

@time fp_arr = map(product(A_range, A_range, A_range, A_range)) do (Aee, Aei, Aie, Aii)
    mods = (α=(0.4, 0.7), Aie=Aie, Aei=Aei, Aee=Aee, Aii=Aii, firing_θI=0.2, blocking_θI=0.5, save_idxs=nothing, save_on=true, saveat=0.1) 
    sim = prototype(; mods...)

    calculate_fixedpoints(sim.model);
end

fp_arr

end;

@show "done"