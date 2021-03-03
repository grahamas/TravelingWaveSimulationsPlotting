using DrWatson
using TravelingWaveSimulations, WilsonCowanModel, 
    TravelingWaveSimulationsPlotting, 
    Simulation73Plotting, Simulation73
using TravelingWaveSimulationsPlotting: _collapse_to_axes
using Dates
using Makie
using AxisIndices

include(projectdir("_drafts/grouped_bar_plot.jl"))

if !@isdefined(refresh_sweep_arrs)
    refresh_sweep_arrs = false
end

if !@isdefined(blocking_fp_arr) || !@isdefined(monotonic_fp_arr) || refresh_sweep_arrs
    A_range = 0.1:0.1:1.5
    sweeping_mods = (Aee=A_range, Aei=A_range, Aie=A_range, Aii=A_range)
    static_mods = (
        α=(0.4, 0.7), 
        firing_θI=0.2, blocking_θI=0.5, 
        n_lattice = 2,
        save_idxs=nothing, save_on=true, saveat=0.1
    ) 
    file, filename = produce_or_load(datadir(), []; 
        prefix = "blocking_fp_arr",
        force = refresh_sweep_arrs
    ) do c
        blocking_fp_arr = sweep_calculate_fixedpoints(
            "full_dynamics_blocking", 
            static_mods,
            sweeping_mods,
            ; 
            dx = 0.01
        )
        return @dict(blocking_fp_arr)
    end
    @unpack blocking_fp_arr = file
    file, filename = produce_or_load(datadir(), []; 
        prefix = "monotonic_fp_arr",
        force = refresh_sweep_arrs
    ) do c
        monotonic_fp_arr = sweep_calculate_fixedpoints(
            "full_dynamics_monotonic", 
            static_mods,
            sweeping_mods,
            ; 
            dx = 0.01
        )
        return @dict(monotonic_fp_arr)
    end
    @unpack monotonic_fp_arr = file
end

refresh_sweep_arrs = false

fig = let example_name = "seizures_fig",
    session_id = "$(Dates.now())",
    example_dir = mkpath(plotsdir("$(example_name)_$(session_id)")),
    stim_strengths = [0.0, 0.02, 0.1],
    figure_resolution = (3000, 1600),
    colorbar_width = 25,
    fp_arrs = (monotonic = monotonic_fp_arr,
        blocking = blocking_fp_arr);

# A_mods = (
#     Aee=1., Aei=0.8,
#     Aie=0.8, Aii=0.2
# )
A_mods = (
    Aee=1.5, Aei=1.5,
    Aie=1.5, Aii=1.5
)
mods = (α=(0.4, 0.7),  
    firing_θI=0.2, θI=0.2, blocking_θI=0.5, 
    save_idxs=nothing, save_on=true, saveat=0.1,
    x_lattice=256., n_lattice=512,
    stim_radius = 1.0,
    stop_time = 4,
    A_mods...
) 
simple_theme = Theme(
    linewidth = 20.0,
    Axis = (
        backgroundcolor = :white,
        leftspinevisible = true,
        rightspinevisible = false,
        bottomspinevisible = true,
        topspinevisible = false,
        xgridcolor = :white,
        ygridcolor = :white,
    )
)

with_theme(simple_theme) do 
    fig = Figure(resolution = figure_resolution)
    for (i_nonl_type, (nonl_sym, fp_arr)) in enumerate(pairs(fp_arrs))
        prototype = get_prototype("full_dynamics_$(string(nonl_sym))");
        sim = prototype(; mods...)

        params = get_nullcline_params(sim)
        # test that the FP counter isn't unstable near this model
        # @assert reduce(==, length.(calculate_fixedpoints.(Ref(params), [0.01, 0.001])))
        # make sure what's saved matches our mods
        calculated_fp = calculate_fixedpoints(params, 0.01)
        saved_fp = getindex(fp_arr; A_mods...)
        # @assert calculated_fp == saved_fp

        fig[i_nonl_type, 1] = plot_nullclines!(fig, params, 0.01; 
            mark_fp = true)

        collapsed_axes = (:Aie, :Aii)
        # min E/I
        ei_differences = map(fp_arr) do fps
            map(fps) do (e_val, i_val)
                e_val - i_val
            end::Vector{Float64}
        end
        @warn "Starting min"
        fig[i_nonl_type, 2] = axisarray_heatmap!(fig, 
            dropdims(minimum(minimum.(ei_differences), dims=collapsed_axes), dims=collapsed_axes), 
            colorbar_width=colorbar_width,
            colorbar_label = "minimum E-I"
        )
        @warn "Starting max"
        fig[i_nonl_type, 3] = axisarray_heatmap!(fig, 
            dropdims(maximum(maximum.(ei_differences), dims=collapsed_axes), dims=collapsed_axes), 
            colorbar_width=colorbar_width,
            colorbar_label = "maximum E-I"
        )
        @warn "Done."

        
    end # for nonl_sym
    fig
end # do with_theme

end # let
