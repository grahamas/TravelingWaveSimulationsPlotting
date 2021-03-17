using DrWatson
using TravelingWaveSimulations, WilsonCowanModel, 
    TravelingWaveSimulationsPlotting, 
    Simulation73Plotting, Simulation73
using TravelingWaveSimulationsPlotting: _collapse_to_axes, calculate_fixedpoints
using Dates
using Makie

let example_name = "he_mono5_block7",
    session_id = "$(Dates.now())",
    example_dir = mkpath(plotsdir("$(example_name)_$(session_id)")),
    figure_resolution = (1200, 700),
    colorbar_width = 25
;

mods = (α=(0.4, 0.7), 
    Aie=0.81, Aei=0.8, 
    firing_θI=0.2, θI=0.2, blocking_θI=0.5, 
    save_idxs=nothing, save_on=true, saveat=0.1,
    x_lattice=256., n_lattice=512,
    stim_radius = 1.0,
    stop_time = 4
) 
simple_theme = Theme(
    fontsize=20,
    linewidth=1.5,
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
        params = get_nullcline_params(sim)

        fig[1, i_nonl_type] = plot_nullclines!(fig, params; linewidth=2.5)
    end # for nonl

    save(
        joinpath(example_dir, 
            "$(example_name).png"
        ), 
        fig
    )

end #with simple_theme

end #let