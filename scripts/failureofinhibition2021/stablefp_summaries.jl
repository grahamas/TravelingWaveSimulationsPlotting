using DrWatson
using TravelingWaveSimulations, WilsonCowanModel, 
    TravelingWaveSimulationsPlotting, 
    Simulation73Plotting, Simulation73
using TravelingWaveSimulationsPlotting: _collapse_to_axes
using Dates
using Makie
#using GLMakie; ext_2d = "png"; GLMakie.activate!()
using CairoMakie; ext_2d = "svg"; CairoMakie.activate!()
using AxisIndices
using AlgebraOfGraphics, ColorSchemes
using DataFrames

# loads failing_fp_arr and monotonic_fp_arr
sub_A_sweep_lower_bound = 0.5
sub_A_sweep_upper_bound = 1.5
sub_A_range = sub_A_sweep_lower_bound..sub_A_sweep_upper_bound
include(scriptsdir("load/he_sweep_arrs.jl"))




let nt_map(fn::Function, nt::NamedTuple{NS}) where NS = NamedTuple{NS}(map(fn, values(nt))),
    fp_arr_nt = (
        fire_fail=blocking_fp_arr[Aee=sub_A_range,Aei=sub_A_range,Aie=sub_A_range,Aii=sub_A_range],
        fire=monotonic_fp_arr[Aee=sub_A_range,Aei=sub_A_range,Aie=sub_A_range,Aii=sub_A_range]
    ),
    nonl_types = keys(fp_arr_nt),
    fp_count_arr_nt = nt_map(x -> length.(x), fp_arr_nt),
    prototype_name_nt = (
        fire="full_dynamics_monotonic",
        fire_fail="full_dynamics_blocking"
    ),
    E_bounds = [0.05, 0.71],
    SI_bounds = [0.05, 0.71],
    arrows_step=0.05,
    smods = (
        α=(0.4, 0.7), 
        firing_θI=0.2, blocking_θI=0.5, 
        n_lattice = 2,
        save_idxs=nothing, save_on=false, saveat=0.1
    ),
    session_name = "stablefp_summaries_$(ext_2d)",
    session_id = "$(Dates.now())",
    axis = (width = 800, height = 800),
    abbrev_count_label(x) = begin
        if x >= 1000
            try
                "$(Int(x / 1000))K"
            catch
                "$(x / 1000)K"
            end
        else
            "$(Int(x))"
        end
    end,
    bar_theme = Theme(
        fontsize=56,
        strokewidth= 5.,
        Axis = (
            backgroundcolor = RGBA(1.,1.,1.,0.),
            leftspinevisible = true,
            rightspinevisible = false,
            bottomspinevisible = true,
            topspinevisible = false,
            xgridcolor = RGBA(1.,1.,1.,0.),
            ygridcolor = RGBA(1.,1.,1.,0.),
            strokewidth= 5.,
            ytickformat = xs -> abbrev_count_label.(xs)
        )
    ),
    nullcline_theme = Theme(
        fontsize=48,
        Axis = (
            backgroundcolor = RGBA(1.,1.,1.,0.),
            leftspinevisible = true,
            rightspinevisible = false,
            bottomspinevisible = true,
            topspinevisible = false,
            strokewidth=2.,
            xgridcolor = RGBA(1.,1.,1.,0.),
            ygridcolor = RGBA(1.,1.,1.,0.)
        ),
        Lines = (
            linewidth=4.0,
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
    )  
;
plots_subdir = "$(session_name)_$(session_id)"
mkpath(plotsdir(plots_subdir))

# # get only stable fixedpoints
# stable_fp_arr_nt = NamedTuple{nonl_types}(
#     map(nonl_types) do nonl_type
#         filter_stable_fps(
#             prototype_name_nt[nonl_type],
#             smods,
#             fp_arr_nt[nonl_type]
#         )
#     end
# )

# calculate seizure index of stable fixed points
stablefp_SI_nt = NamedTuple{nonl_types}(
    map(nonl_types) do nonl_type
        prototype = get_prototype(prototype_name_nt[nonl_type])
        fp_arr = fp_arr_nt[nonl_type]
        stablefp_SI = Array{Union{Float64,Missing}}(missing, size(fp_arr)..., 7)
        for (mx_idx, (nt, fps)) ∈ zip(CartesianIndices(fp_arr), TravelingWaveSimulationsPlotting.enumerate_nt(fp_arr))
            params = get_nullcline_params(prototype(; smods..., nt...))
            for fp_idx ∈ 1:length(fps)
                fp = fps[fp_idx]
                if TravelingWaveSimulationsPlotting.fixedpoint_is_stable(params, fp)
                    stablefp_SI[mx_idx, fp_idx] = seizure_index(fp)
                end
            end
        end
        stablefp_SI
    end
)
stablefp_E_nt = NamedTuple{nonl_types}(
    map(nonl_types) do nonl_type
        prototype = get_prototype(prototype_name_nt[nonl_type])
        fp_arr = fp_arr_nt[nonl_type]
        stablefp_E = Array{Union{Float64,Missing}}(missing, size(fp_arr)..., 7)
        for (mx_idx, (nt, fps)) ∈ zip(CartesianIndices(fp_arr), TravelingWaveSimulationsPlotting.enumerate_nt(fp_arr))
            params = get_nullcline_params(prototype(; smods..., nt...))
            for fp_idx ∈ 1:length(fps)
                fp = fps[fp_idx]
                if TravelingWaveSimulationsPlotting.fixedpoint_is_stable(params, fp)
                    stablefp_E[mx_idx, fp_idx] = first(fp)
                end
            end
        end
        stablefp_E
    end
)

n_obs = length(stablefp_E_nt[nonl_types[1]])
unrolled_E_df = DataFrame(
    nonl_type=repeat(nonl_types |> collect, inner=(n_obs,)),
    E = vcat([stablefp_E_nt[nonl][:] for nonl in nonl_types]...)
)

unrolled_SI_df = DataFrame(
    nonl_type=repeat(nonl_types |> collect, inner=(n_obs,)),
    SI = vcat([stablefp_SI_nt[nonl][:] for nonl in nonl_types]...)
)
dropmissing!(unrolled_E_df); dropmissing!(unrolled_SI_df)
unrolled_E_df.condition = map(unrolled_E_df.E) do E
    if E .< E_bounds[1]
        "near zero"
    elseif  E_bounds[1] .< E .< E_bounds[2]
        "not extreme"#"$(E_bounds[1]) < E < $(E_bounds[2])"
    else
        "∼$(E_bounds[2])"
    end
end

unrolled_SI_df.condition = map(unrolled_SI_df.SI) do SI
    if SI .< SI_bounds[1]
        "near zero"
    elseif  SI_bounds[1] .< SI .< SI_bounds[2]
        "not extreme"#"$(SI_bounds[1]) < SI < $(SI_bounds[2])"
    else
        "∼$(SI_bounds[2])"
    end
end

with_theme(bar_theme) do
    SI_condition_frequency = data(unrolled_SI_df) * frequency() * mapping(:nonl_type)
    SI_stacks = SI_condition_frequency * mapping(color=:condition, stack=:condition)
    axis = merge(axis, (title="seizure index",))
    fig = draw(SI_stacks; axis)
    fig.figure.current_axis.x.xlabel = "inh. nonlinearity"
    tightylimits!(fig)
    save(plotsdir(plots_subdir, "stablefp_SI_distribution.$(ext_2d)"), fig)

    E_condition_frequency = data(unrolled_E_df) * frequency() * mapping(:nonl_type)
    E_stacks = E_condition_frequency * mapping(color=:condition, stack=:condition)
    axis = merge(axis, (title="excitatory activity",))
    fig = draw(E_stacks; axis)
    fig.figure.current_axis.x.xlabel = "inh. nonlinearity"
    tightylimits!(fig)
    save(plotsdir(plots_subdir, "stablefp_E_distribution.$(ext_2d)"), fig)

    fig

end # bar_theme

# with_theme(nullcline_theme) do
#     fig = Figure()
#     fig[1,1] = ax = Makie.Axis(fig)
#     xs = ones(Int, length(unrolled_SI_df.nonl_type))
#     xs[unrolled_SI_df.nonl_type .== nonl_types[2]] .= 2
#     violin!(ax, xs, unrolled_SI_df.SI, width=0.5)
#     ax.xticks = [1,2]
#     ax.xtickformat = xs -> ["monotonic", "failing"]
#     fig
# end


end # let