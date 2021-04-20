using WilsonCowanModel: wcm_du_defn, wcm_dv_defn, WCMParams, AbstractNullclineParams, get_nullcline_params
using Makie: @lift


function lifted_wcm_param(;
    Aee=1., Aie=1., Aei=1.5, Aii=0.25,
    θef=0.125, θif=0.4, θib=7., τ=0.4, β=50.,
    decaye=1., decayi=1.)
    @lift WCMParams(;
            Aee=$(Node(Aee)), Aie=$(Node(Aie)), Aei=$(Node(Aei)), Aii=$(Node(Aii)),
            θef=$(Node(θef)), θif=$(Node(θif)), θib=$(Node(θib)), τ=$(Node(τ)), β=$(Node(β)),
            decaye=$(Node(decaye)), decayi=$(Node(decayi))
        )
end

using OffsetArrays

function plot_nullclines!(fig::Figure, p::Union{AbstractWCMDepNullclineParams, AbstractWCMNullclineParams}, axis_length::Integer=100; kwargs...)
    plot_nullclines!(fig, p, [range(0., 1., length=axis_length), range(0., 1., length=axis_length)]; kwargs...)
end
function plot_nullclines!(fig::Figure, p::Union{AbstractWCMDepNullclineParams, AbstractWCMNullclineParams}, nullcline_axes::AbstractVector{<:AbstractVector};
        xlabel="E", ylabel="I", 
        mark_fp=true, arrows_step=nothing,
        linewidth=1)
    ax = MakieLayout.Axis(fig, aspect=DataAspect(),
        xlabel=xlabel, ylabel=ylabel
    )

    dus = calculate_field(wcm_du_defn, nullcline_axes, p)
    dvs = calculate_field(wcm_dv_defn, nullcline_axes, p)
    u_nullclines = lines(contour(nullcline_axes..., dus, 0.))
    v_nullclines = lines(contour(nullcline_axes..., dvs, 0.))

    for line in u_nullclines
        xs, ys = coordinates(line)
        Makie.lines!(ax, xs, ys, color=:blue, linewidth=linewidth)
    end

    for line in v_nullclines
        xs, ys = coordinates(line)
        Makie.lines!(ax, xs, ys, color=:red, linestyle=:dash, linewidth=linewidth)
    end

    if mark_fp
        fixedpoints = calculate_fixedpoints(p, length(nullcline_axes[1]))
        stability = fixedpoint_stability.(Ref(p), fixedpoints)
        possible_stability_markers = OffsetArray([:xcross, :star, :circle], -1:1)
        stability_marker = getindex.(Ref(possible_stability_markers), stability)
        scatter!(ax, Point2f0.(fixedpoints), marker=stability_marker, markersize=15)
    end

    if arrows_step !== nothing
        arrow_axis = 0.:arrows_step:1. |> collect
        arrow_axes = [arrow_axis, arrow_axis]
        dus = calculate_field(wcm_du_defn, arrow_axes, p)
        dvs = calculate_field(wcm_dv_defn, arrow_axes, p)
        arrows!(ax, arrow_axis, arrow_axis, dus, dvs; arrowsize = 0.03, normalize=true, lengthscale=0.02f0)
    end

    xlims!(ax, 0., 1.)
    ylims!(ax, 0., 1.)

    return ax
end




function calculate_fixedpoints(mdb_path::AbstractString, new_mods::NamedTuple{NAMES}, dx::T=0.01) where {NAMES, T<:Number}
    prototype_name, saved_mods = read_params_from_data_path(mdb_path)
    @assert NAMES ⊆ keys(saved_mods)
    combined_mods = merge(saved_mods, Dict(pairs(new_mods))) # overwrites saved_mods
    
    variable_mod_names, variable_mod_values = zip([(name, values) for (name,values) in pairs(combined_mods) if length(values)!= 1]...)
    fixedpoints_arr = NamedAxisArray{variable_mod_names}(
        Array{Vector{SVector{2,T}},length(variable_mod_values)}(
            undef, length.(variable_mod_values)...
        ), 
        variable_mod_values...
    )

    prototype = get_prototype(prototype_name)

    us = 0.:dx:1.
    vs = copy(us)
    dus = Array{T,2}(undef, length(us), length(vs))

    map(product(variable_mod_values...)) do mod_values
        fixing_mods = NamedTuple{variable_mod_names}(mod_values)
        # other_opts = Dict() makes sure it saves all frames
        mods = (combined_mods..., fixing_mods..., save_idxs=nothing, save_on=true)
        #write_modifications!(plots_subdir, mods, unique_id)
        model = prototype(; mods...)
        setindex!(fixedpoints_arr, calculate_fixedpoints!(dus, (us, vs), model);
            fixing_mods...)
    end
    return fixedpoints_arr, combined_mods, prototype
end

using Roots
function Roots.find_zero(fn::Function, (p1,p2)::Tuple{PT,PT}, args...; kwargs...) where {N,T, PT<:SVector{N,T}}
    interp(x) = x .* p1 .+ (1-x) .* p2
    interp_zero = find_zero(x -> fn(interp(x)), (0, 1), args...; kwargs...)
    return interp(interp_zero)
end

function edge_push_intersections!(intersections, prealloc_subgrid,
        edge_vertex, interior_vertex, nullcline_params, 
        supergrid_dx, subgrid_n, subgrid_x, fp_zero_atol)
    """ On the edge, don't insist on dV crossing zero; might just get close."""
    dist_between_points = _dist(edge_vertex, interior_vertex)
    subgrid_center_dists = 0.:supergrid_dx:dist_between_points # step at same res as supergrid
    d_vertex = interior_vertex .- edge_vertex 
    d_vertex /= norm(d_vertex)
    @assert isapprox(edge_vertex + (dist_between_points * d_vertex), interior_vertex, atol=sqrt(eps())) "$(edge_vertex + (dist_between_points * d_vertex)) != $(interior_vertex)"
    # FIXME delete assert ^ after this works
    edge_subgrid_push_intersections!(intersections, prealloc_subgrid,
        edge_vertex, nullcline_params, subgrid_n, subgrid_x, fp_zero_atol
    )
    for subgrid_center_dist ∈ subgrid_center_dists[begin:end]
        subgrid_center = edge_vertex + (subgrid_center_dist * d_vertex)
        interior_subgrid_push_intersections!(intersections,
            prealloc_subgrid, subgrid_center,
            nullcline_params, 
            subgrid_n, subgrid_x, fp_zero_atol
        )
    end
end

function edge_subgrid_push_intersections!(intersections, prealloc_subgrid, point::SVector{2,Float64}, params, subgrid_n::Int, subgrid_x, fp_zero_atol)
    # calculate dv into subgrid
    subgrid_us, subgrid_vs = calculate_dv_subgrid!(prealloc_subgrid, point, params, subgrid_n, subgrid_x)
    # assume there can only be one fixed point in the subgrid
    # first check for interior fixed point
    pre_count = length(intersections)
    _interior_subgrid_push_intersections!(intersections, subgrid_us, subgrid_vs, prealloc_subgrid, params, subgrid_x/subgrid_n, fp_zero_atol)
    if length(intersections) == pre_count # no interior intersection
        edge_min_val = Inf
        edge_min_vertex = nothing
        # the first entry of subgrid_us and/or subgrid_vs will define the edge
        @assert subgrid_us[end] != 0 && subgrid_us[begin] != 1 && subgrid_vs[end] != 0 && subgrid_vs[begin] != 1 "subgrid_us: $subgrid_us \n subgrid_vs: $subgrid_vs"
        # FIXME should really check du
        if subgrid_us[begin] == 0 || subgrid_us[begin] == 1
            (min_val, min_idx) = findmin(abs.(prealloc_subgrid[begin,:]) .+ abs.(wcm_du_defn.(subgrid_us[begin], subgrid_vs, Ref(params))))
            if min_val < fp_zero_atol
                edge_min_val = min_val
                edge_min_vertex = SVector{2,Float64}(subgrid_us[begin], subgrid_vs[min_idx])
            end
        end
        if subgrid_us[end] == 0 || subgrid_us[end] == 1
        (min_val, min_idx) = findmin(abs.(prealloc_subgrid[begin,:]) .+ abs.(wcm_du_defn.(subgrid_us[end], subgrid_vs, Ref(params))))
            if min_val < fp_zero_atol
                edge_min_val = min_val
                edge_min_vertex = SVector{2,Float64}(subgrid_us[end], subgrid_vs[min_idx])
            end
        end
        if subgrid_vs[begin] == 0 || subgrid_vs[begin] == 1
            (min_val, min_idx) = findmin(abs.(prealloc_subgrid[begin,:]) .+ abs.(wcm_du_defn.(subgrid_us, subgrid_vs[begin], Ref(params))))
            if min_val < edge_min_val && min_val < fp_zero_atol
                edge_min_vertex = SVector{2,Float64}(subgrid_us[min_idx], subgrid_vs[begin])
            end
        end
        if subgrid_vs[end] == 1 || subgrid_vs[end] == 0
        (min_val, min_idx) = findmin(abs.(prealloc_subgrid[begin,:]) .+ abs.(wcm_du_defn.(subgrid_us, subgrid_vs[end], Ref(params))))
            if min_val < edge_min_val && min_val < fp_zero_atol
                edge_min_vertex = SVector{2,Float64}(subgrid_us[min_idx], subgrid_vs[end])
            end
        end
        if edge_min_vertex !== nothing
            push_intersection!(intersections, edge_min_vertex, params, subgrid_x/subgrid_n, fp_zero_atol)
        end
    end
end
  
function nt_index(arr::NamedAxisArray{nt_names}, idx) where nt_names
    nt_vals = collect(product(axes_keys(arr)...))[idx]
    NamedTuple{nt_names}(nt_vals)
end
