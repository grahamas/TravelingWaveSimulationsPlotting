using DrWatson
using TravelingWaveSimulations, WilsonCowanModel, 
    TravelingWaveSimulationsPlotting, 
    Simulation73Plotting, Simulation73
using Dates
using Makie

# if !@isdefined(blocking_fp_arr)
#     @warn "Calculating blocking_fp_arr..."
#     include(scriptsdir("he_blocking_fp_sweep.jl"))
#     @warn "done."
# end

# if !@isdefined(mono_fp_arr)
#     @warn "Calculating monotonic_fp_arr..."
#     include(scriptsdir("he_monotonic_fp_sweep.jl"))
#     @warn "done."
# end

let example_name = "he_5vs7fp",
    session_id = "$(Dates.now())",
    example_dir = mkpath(plotsdir("$(example_name)_$(session_id)")),
    stim_strengths = [0.0, 0.02, 0.1],
    figure_resolution = (3000, 1600);

mods = (α=(0.4, 0.7), 
    Aie=0.81, Aei=0.8, 
    firing_θI=0.2, θI=0.2, blocking_θI=0.5, 
    save_idxs=nothing, save_on=true, saveat=0.1,
    x_lattice=256., n_lattice=512,
    stim_radius = 1.0,
    stop_time = 4
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
nonl_type_strings = ["monotonic", "blocking"]

with_theme(simple_theme) do 
    fig = Figure(resolution = figure_resolution)
    for (i_nonl_type, nonl_type_string) ∈ enumerate(nonl_type_strings)
        prototype = get_prototype("full_dynamics_$(nonl_type_string)");
        sim = prototype(; mods...)

        # test that the FP counter isn't unstable near this model
        @assert reduce(==, calculate_fixedpoints.(Ref(sim.model), [0.01, 0.001]))
        params = get_nullcline_params(sim)
        fig[i_nonl_type, 1] = plot_nullclines!(fig, params, 0.01)
        execs_by_stim_strength = [
            execute(prototype(; mods..., stim_strength=stim_strength)) 
            for stim_strength ∈ stim_strengths
        ]
        fig[i_nonl_type, 1 .+ (1:length(stim_strengths))] = 
            isolimit_exec_heatmaps!(
                fig, execs_by_stim_strength
                ; titles = ["S = $stim_str" for stim_str ∈ stim_strengths]
        );
    end # for nonl_type
    save(
        joinpath(example_dir, 
            "$(example_name).png"
        ), 
        fig
    )
end # with_theme
end # let
