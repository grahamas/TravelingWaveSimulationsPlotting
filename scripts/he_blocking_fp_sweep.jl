using TravelingWaveSimulationsPlotting, TravelingWaveSimulations
using WilsonCowanModel, Simulation73
using Dates
using Contour
using AxisIndices, IterTools
using Makie

@show "running"


blocking_fp_arr = let nonl_type = "blocking",
    session_name = "he_$(nonl_type)_sweep_fp",
    session_id = "$(Dates.now())";

# 7fp values:
# Aee = 1.
# Aei = 0.8
# Aie = 1.
# Aii = 0.25
A_range = 0.1:0.1:1.5
@warn "Checking $(length(A_range)^4) parameterizations..."

prototype = get_prototype("full_dynamics_$(nonl_type)")

@time fp_arr = map(product(A_range, A_range, A_range, A_range)) do (Aee, Aei, Aie, Aii) # ATTN gives ordering of names in NamedAxisArray defn below
    mods = (α=(0.4, 0.7), Aie=Aie, Aei=Aei, Aee=Aee, Aii=Aii, firing_θI=0.2, blocking_θI=0.5, save_idxs=nothing, save_on=true, saveat=0.1) 
    sim = prototype(; mods...)

    calculate_fixedpoints(sim.model);
end

fp_count_arr = NamedAxisArray{(:Aee, :Aei, :Aie, :Aii)}(length.(fp_arr), A_range, A_range, A_range, A_range) # ATTN order comes from map arg destructuring above

n_fps = 0:7
fp_count = [count(fp_count_arr .== x) for x in n_fps]
max_log = log10.(fp_count) |> maximum |> mx -> ceil(Int, mx)

sc, ly = layoutscene(); ly[1,1] = ax = MakieLayout.Axis(sc);
barplot!(ax, n_fps, log10.(fp_count))
tightlimits!(ax)
xlims!(ax, n_fps[begin]-0.5,n_fps[end]+0.5)
ylims!(ax, 0, max_log)
ax.xticks=n_fps
ax.yticks=0:max_log
ax.ytickformat = xs -> [x > 0 ? "10^$(Int(x))" : "$x" for x in xs]
ax.title = "Varying connectivity of $(nonl_type) WCM"
ax.xlabel = "# fixed points"
ax.ylabel = "# models"
hidespines!(ax, :t, :r)
hidedecorations!(ax, label=false, ticklabels=false, ticks=false)
mkpath(plotsdir("$(session_name)_$(session_id)"))
save(plotsdir("$(session_name)_$(session_id)", "he_$(nonl_type)_fp.png"), sc)
 
fp_arr

end;

@show "done"