using DrWatson
using TravelingWaveSimulations, WilsonCowanModel, 
    TravelingWaveSimulationsPlotting, 
    Simulation73Plotting, Simulation73
using TravelingWaveSimulationsPlotting: _collapse_to_axes
using Dates
using Makie
using AxisIndices

# loads blocking_fp_arr and monotonic_fp_arr

let blocking_prototype_name = "full_dynamics_blocking",
    monotonic_prototype_name = "full_dynamics_monotonic",
    smods = (
        α=(0.4, 0.7), 
        firing_θI=0.2, blocking_θI=0.5, 
        n_lattice = 2,
        save_idxs=nothing, save_on=false, saveat=0.1
    ),
    session_name = "pharmacology_phase_spaces",
    session_id = "$(Dates.now())",
    figure_resolution=(800,800),
    simple_theme = Theme(
        linewidth = 20.0,
        fontsize=24,
        Axis = (
            backgroundcolor = :white,
            leftspinevisible = true,
            rightspinevisible = false,
            bottomspinevisible = true,
            topspinevisible = false,
            xgridcolor = :white,
            ygridcolor = :white,
        )
    ),
    plots_subdir = plotsdir("$(session_name)_$(session_id)")
;
mkpath(plots_subdir)

vulnerable_3sfp = (Aee = 1.2, Aei = 1.0, Aie = 1.4, Aii = 0.8)
low_Aie_3sfp = (Aee = 1.2, Aei = 1.0, Aie = 0.9, Aii = 0.8)
low_Aee_3sfp = (Aee = 0.9, Aei = 1.0, Aie = 1.4, Aii = 0.8)
low_excitation_3sfp = (Aee = 0.9, Aei = 1.0, Aie = 0.9, Aii = 0.8)
high_Aii_3sfp = (Aee = 1.2, Aei = 1.0, Aie = 1.4, Aii = 1.3)
high_Aei_3sfp = (Aee = 1.2, Aei = 1.4, Aie = 1.4, Aii = 0.8)
high_inhibition_3sfp = (Aee = 1.2, Aei = 1.4, Aie = 1.4, Aii = 1.3)

for (sym, prototype_name) ∈ [(:blocking, blocking_prototype_name), (:monotonic, monotonic_prototype_name)]
    prototype = get_prototype(prototype_name)

    plot_and_save(plot_nullclines!,
        get_nullcline_params(prototype(; smods..., vulnerable_3sfp...));
        plot_name="nullclines_$(sym)_vulnerable_3sfp.png",
        plots_subdir=plots_subdir,
        arrows_step=0.05,
        title="$(vulnerable_3sfp)"
    )
    
    plot_and_save(plot_nullclines!,
        get_nullcline_params(prototype(; smods..., low_Aie_3sfp...));
        plot_name="nullclines_$(sym)_low_Aie_3sfp.png",
        plots_subdir=plots_subdir,
        arrows_step=0.05,
        title="$(low_Aie_3sfp)"
    )

    plot_and_save(plot_nullclines!,
        get_nullcline_params(prototype(; smods..., low_Aee_3sfp...));
        plot_name="nullclines_$(sym)_low_Aee_3sfp.png",
        plots_subdir=plots_subdir,
        arrows_step=0.05,
        title="$(low_Aee_3sfp)"
    )

    plot_and_save(plot_nullclines!,
        get_nullcline_params(prototype(; smods..., low_excitation_3sfp...));
        plot_name="nullclines_$(sym)_low_excitation_3sfp.png",
        plots_subdir=plots_subdir,
        arrows_step=0.05,
        title="$(low_excitation_3sfp)"
    )

    plot_and_save(plot_nullclines!,
        get_nullcline_params(prototype(; smods..., high_Aii_3sfp...));
        plot_name="nullclines_$(sym)_high_Aii_3sfp.png",
        plots_subdir=plots_subdir,
        arrows_step=0.05,
        title="$(high_Aii_3sfp)"
    )

    plot_and_save(plot_nullclines!,
        get_nullcline_params(prototype(; smods..., high_Aei_3sfp...));
        plot_name="nullclines_$(sym)_high_Aei_3sfp.png",
        plots_subdir=plots_subdir,
        arrows_step=0.05,
        title="$(high_Aei_3sfp)"
    )

    plot_and_save(plot_nullclines!,
        get_nullcline_params(prototype(; smods..., high_inhibition_3sfp...));
        plot_name="nullclines_$(sym)_high_inhibition_3sfp.png",
        plots_subdir=plots_subdir,
        arrows_step=0.05,
        title="$(high_inhibition_3sfp)"
    )
    
    
end



end