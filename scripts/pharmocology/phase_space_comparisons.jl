using DrWatson
using TravelingWaveSimulations, WilsonCowanModel, 
    TravelingWaveSimulationsPlotting, 
    Simulation73Plotting, Simulation73
using TravelingWaveSimulationsPlotting: _collapse_to_axes
using Dates
#using GLMakie; ext_2d = "png"; GLMakie.activate!(); @warn "Using OpenGL"
using CairoMakie; ext_2d = "svg"; CairoMakie.activate!(); @warn "Using Cairo"
using ColorSchemes

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
    session_name = "pharmacology_phase_spaces_arrows_$(ext_2d)",
    session_id = "$(Dates.now())",
    figure_resolution=(800,800),
    simple_theme = Theme(
        Lines = (
            linewidth=4.0,
        ),
        fontsize=30,
        Axis = (
            backgroundcolor = :white,
            leftspinevisible = true,
            rightspinevisible = false,
            bottomspinevisible = true,
            topspinevisible = false,
            xgridcolor = :white,
            ygridcolor = :white,
        ),
        Arrows = (
            arrowsize=10, lengthscale=0.017,
            linewidth=2,
            arrowcolor=:black, linecolor=:black,
            colormap=ColorSchemes.Greys_5,
            normalize=false
        ),
        Scatter = (
            markersize=27,
            strokewidth=1
        )
    ),
    plots_subdir = plotsdir("$(session_name)_$(session_id)"),
    arrows_step=0.05
;
mkpath(plots_subdir)

vulnerable_3sfp = (Aee = 1.2, Aei = 1.0, Aie = 1.4, Aii = 0.8)
low_Aie_3sfp = (Aee = 1.2, Aei = 1.0, Aie = 0.9, Aii = 0.8)
low_Aee_3sfp = (Aee = 0.9, Aei = 1.0, Aie = 1.4, Aii = 0.8)
low_excitation_3sfp = (Aee = 0.9, Aei = 1.0, Aie = 0.9, Aii = 0.8)
high_Aii_3sfp = (Aee = 1.2, Aei = 1.0, Aie = 1.4, Aii = 1.3)
high_Aei_3sfp = (Aee = 1.2, Aei = 1.4, Aie = 1.4, Aii = 0.8)
high_inhibition_3sfp = (Aee = 1.2, Aei = 1.4, Aie = 1.4, Aii = 1.3)

θef = 0.125
θif = 0.2
θib = 0.5
Δθ = 0.2

with_theme(simple_theme) do
for (sym, prototype_name) ∈ [(:blocking, blocking_prototype_name), (:monotonic, monotonic_prototype_name)]
    prototype = get_prototype(prototype_name)

    plot_and_save(plot_nullclines!,
        get_nullcline_params(prototype(; smods..., vulnerable_3sfp...));
        plot_name="nullclines_$(sym)_vulnerable_3sfp_$(vulnerable_3sfp).$(ext_2d)",
        plots_subdir=plots_subdir,
        arrows_step=arrows_step
    )
    
    plot_and_save(plot_nullclines!,
        get_nullcline_params(prototype(; smods..., low_Aie_3sfp...));
        plot_name="nullclines_$(sym)_low_Aie_3sfp_$(low_Aie_3sfp).$(ext_2d)",
        plots_subdir=plots_subdir,
        arrows_step=arrows_step
    )

    # plot_and_save(plot_nullclines!,
    #     get_nullcline_params(prototype(; smods..., low_Aee_3sfp...));
    #     plot_name="nullclines_$(sym)_low_Aee_3sfp_$(low_Aee_3sfp).$(ext_2d)",
    #     plots_subdir=plots_subdir,
    #     arrows_step=arrows_step,
    #     #title="ΔAee = $(round(low_Aee_3sfp[:Aee] - vulnerable_3sfp[:Aee], digits=1))"
    # )

    # plot_and_save(plot_nullclines!,
    #     get_nullcline_params(prototype(; smods..., low_excitation_3sfp...));
    #     plot_name="nullclines_$(sym)_low_excitation_3sfp_$(low_excitation_3sfp).$(ext_2d)",
    #     plots_subdir=plots_subdir,
    #     arrows_step=arrows_step,
    #     #title="ΔAee = $(round(low_excitation_3sfp[:Aee] - vulnerable_3sfp[:Aee], digits=1)); ΔAie = $(round(low_excitation_3sfp[:Aie] - vulnerable_3sfp[:Aie], digits=1))"
    # )

    # plot_and_save(plot_nullclines!,
    #     get_nullcline_params(prototype(; smods..., high_Aii_3sfp...));
    #     plot_name="nullclines_$(sym)_high_Aii_3sfp_$(high_Aii_3sfp).$(ext_2d)",
    #     plots_subdir=plots_subdir,
    #     arrows_step=arrows_step,
    #     #title="ΔAii = $(round(high_Aii_3sfp[:Aii] - vulnerable_3sfp[:Aii], digits=1))"
    # )

    # plot_and_save(plot_nullclines!,
    #     get_nullcline_params(prototype(; smods..., high_Aei_3sfp...));
    #     plot_name="nullclines_$(sym)_high_Aei_3sfp_$(high_Aei_3sfp).$(ext_2d)",
    #     plots_subdir=plots_subdir,
    #     arrows_step=arrows_step,
    #     #title="ΔAei = $(round(high_Aei_3sfp[:Aei] - vulnerable_3sfp[:Aei], digits=1))"
    # )

    # plot_and_save(plot_nullclines!,
    #     get_nullcline_params(prototype(; smods..., high_inhibition_3sfp...));
    #     plot_name="nullclines_$(sym)_high_inhibition_3sfp_$(high_inhibition_3sfp).$(ext_2d)",
    #     plots_subdir=plots_subdir,
    #     arrows_step=arrows_step,
    #     #title="ΔAii = $(round(high_inhibition_3sfp[:Aii] - vulnerable_3sfp[:Aii], digits=1)); ΔAei = $(round(high_inhibition_3sfp[:Aei] - vulnerable_3sfp[:Aei], digits=1))"
    # )

    # plot_and_save(plot_nullclines!,
    #     get_nullcline_params(prototype(; smods..., vulnerable_3sfp..., firing_θI= θif+Δθ, θE=θef+Δθ));
    #     plot_name="nullclines_$(sym)_EI_firing_Δθ_$(Δθ)_3sfp.$(ext_2d)",
    #     plots_subdir=plots_subdir,
    #     arrows_step=arrows_step,
    #     #title="Firing Δθ: $(Δθ)"
    # )

    # Δθ = -0.005
    # plot_and_save(plot_nullclines!,
    #     get_nullcline_params(prototype(; smods..., vulnerable_3sfp..., θE=θef+Δθ));
    #     plot_name="nullclines_$(sym)_E_firing_Δθ_$(Δθ)_3sfp.$(ext_2d)",
    #     plots_subdir=plots_subdir,
    #     arrows_step=arrows_step,
    #     #title="Firing ΔθE: $(Δθ)"
    # )

    # Δθ = -0.025
    # plot_and_save(plot_nullclines!,
    #     get_nullcline_params(prototype(; smods..., vulnerable_3sfp..., θE=θef+Δθ));
    #     plot_name="nullclines_$(sym)_E_firing_Δθ_$(Δθ)_3sfp.$(ext_2d)",
    #     plots_subdir=plots_subdir,
    #     arrows_step=arrows_step,
    #     #title="Firing ΔθE: $(Δθ)"
    # )

    # Δθ = 0.25
    # plot_and_save(plot_nullclines!,
    #     get_nullcline_params(prototype(; smods..., vulnerable_3sfp..., blocking_θI=θib+Δθ, firing_θI=θif+Δθ));
    #     plot_name="nullclines_$(sym)_I_both_Δθ_$(Δθ)_3sfp.$(ext_2d)",
    #     plots_subdir=plots_subdir,
    #     arrows_step=arrows_step,
    #     #title="Both ΔθI: $(Δθ)"
    # )

    # if prototype_name == blocking_prototype_name
    #     Δθ = 0.3
    #     plot_and_save(plot_nullclines!,
    #         get_nullcline_params(prototype(; smods..., vulnerable_3sfp..., blocking_θI=θib+Δθ));
    #         plot_name="nullclines_$(sym)_I_blocking_Δθ_$(Δθ)_3sfp.$(ext_2d)",
    #         plots_subdir=plots_subdir,
    #         arrows_step=arrows_step,
    #         #title="Blocking ΔθI: $(Δθ)"
    #     )
    # end
    

end # for
end # with_theme
end # let