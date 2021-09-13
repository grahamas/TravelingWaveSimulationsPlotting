# using DrWatson
# using TravelingWaveSimulations, WilsonCowanModel,
#     Simulation73Plotting, Simulation73
# using Dates
# using Makie
# using AxisIndices, IterTools

using ProgressMeter

function wcm_sweep_calculate_fixedpoints(prototype_name::String, static_mods, dynamic_mods::NamedTuple{NAMES}; axis_length::Integer=150) where {NAMES}
    if :n_lattice ∉ NAMES
        # @warn "Don't forget to make n_lattice small!"
    end
    prototype = get_prototype(prototype_name)
    sweep_axes = values(dynamic_mods)

    us = range(0., 1., length=axis_length)
    vs = copy(us)
    dus = Array{Float64,2}(undef, length(us), length(vs))

    sweep = product(sweep_axes...)
    progress = Progress(length(sweep); dt=1, desc="Calculating fixed points array...", showspeed=true)

    NamedAxisArray{NAMES}(map() do (sweeping_vals...)
        sweeping_mods = NamedTuple{NAMES}(sweeping_vals...)
        sim = prototype(; static_mods..., sweeping_mods...)
        params = get_nullcline_params(sim.model)
        fps = calculate_fixedpoints!(dus, [us, vs], (wcm_du_defn, wcm_dv_defn), params, ((0.,1.), (0.,1.)))
        @assert length(fps) <= 7
        ProgressMeter.next!(progress)
    end, sweep_axes)
end