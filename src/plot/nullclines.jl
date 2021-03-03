using WilsonCowanModel: wcm_du_defn, wcm_dv_defn, WCMParams, AbstractNullclineParams, get_nullcline_params
using Makie: @lift

using DataStructures: MutableLinkedList, ListNode, length


Base.iterate(l::MutableLinkedList) = l.len == 0 ? nothing : (l.node.next, l.node.next.next)
Base.iterate(l::MutableLinkedList, n::ListNode) = n === l.node ? nothing : (n, n.next)
collect_array(l::MutableLinkedList) = [node.data for node in l]


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

using LoopVectorization

function calculate_field!(field_vals::AbstractMatrix{T}, field_fn, us::AbstractVector{T}, vs::AbstractVector{T}, p) where T
    us = collect(us); vs = collect(vs) # FIXME remove line when VectorizationBase updates to v0.18.3(ish)
    @avx for u_idx ∈ axes(us,1), v_idx ∈ axes(vs,1)
        field_vals[u_idx, v_idx] = field_fn(us[u_idx], vs[v_idx], p)
    end
    field_vals
end

function calculate_field(field_fn::Function, us::AbstractVector{T}, vs::AbstractVector{T}, p) where T
    dus = Array{T,2}(undef, length(us), length(vs))
    calculate_field!(dus, field_fn, us, vs, p)
    return dus
end

function plot_nullclines!(fig::Figure, p::AbstractNullclineParams, dx=0.01; kwargs...)
    plot_nullclines!(fig, p, 0.:dx:1., 0.:dx:1.; kwargs...)
end
function plot_nullclines!(fig::Figure, p::AbstractNullclineParams, us::AbstractVector, vs::AbstractVector;
        xlabel="u", ylabel="v", mark_fp=true)
    ax = MakieLayout.Axis(fig, aspect=DataAspect(),
        xlabel=xlabel, ylabel=ylabel
    )

    # dus = [wcm_du_defn(u, v, p) for u in us, v in vs]
    # dvs = [wcm_dv_defn(u, v, p) for u in us, v in vs]
    dus = calculate_field(wcm_du_defn, us, vs, p)
    dvs = calculate_field(wcm_dv_defn, us, vs, p)
    u_nullclines = lines(contour(us, vs, dus, 0.))
    v_nullclines = lines(contour(us, vs, dvs, 0.))

    for line in u_nullclines
        xs, ys = coordinates(line)
        Makie.lines!(ax, xs, ys, color=:blue)
    end

    for line in v_nullclines
        xs, ys = coordinates(line)
        Makie.lines!(ax, xs, ys, color=:red, linestyle=:dash, linesize=5)
    end

    if mark_fp
        fixedpoints = calculate_fixedpoints(u_nullclines, p, step(us))
        scatter!(ax, Point2f0.(fixedpoints))
    end

    return ax
end

export calculate_fixedpoints!

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
        setindex!(fixedpoints_arr, calculate_fixedpoints!(dus, us, vs, model);
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

calculate_fixedpoints(model::Union{AbstractModel,AbstractSimulation}, args...; kwargs...) = calculate_fixedpoints(get_nullcline_params(model), args...; kwargs...)
calculate_fixedpoints!(du_arr::AbstractArray, us::AbstractVector, vs::AbstractVector, model::Union{AbstractModel{T},AbstractSimulation{T}}, args...; kwargs...) where T = calculate_fixedpoints!(du_arr, us, vs, get_nullcline_params(model), args...; kwargs...)

function calculate_fixedpoints(nullcline_params::AbstractNullclineParams, phase_dx::T=0.01; kwargs...) where {T <: Number}
    us = 0.0:phase_dx:1.0
    vs = copy(us)
    dus = Array{T,2}(undef, length(us), length(vs))

    calculate_fixedpoints!(dus, us, vs, nullcline_params; kwargs...)
end

function calculate_fixedpoints!(dus, us, vs, nullcline_params; kwargs...)
    calculate_field!(dus, wcm_du_defn, 
        us, vs, 
        nullcline_params
    )
    u_nullclines = lines(contour(us, vs, dus, 0.))
    calculate_fixedpoints(u_nullclines, nullcline_params, step(us); kwargs...)
end

function calculate_fixedpoints(u_nullclines::Vector{Curve2{T}},
        nullcline_params::AbstractNullclineParams, 
        supergrid_dx;
        subgrid_side_n=20,
        fp_zero_atol=1e-5
    ) where T
    # only need one variable's nullclines bc all FP must be on one of these nullclines
    subgrid_side_x = 1.2supergrid_dx
    intersections = MutableLinkedList{SVector{2,T}}()
    prealloc_subgrid = Array{T,2}(undef, subgrid_side_n, subgrid_side_n)
    for u_line ∈ u_nullclines
        vertices = u_line.vertices
        @assert length(vertices) ≥ 3
        first_vertex = vertices[begin]
        second_vertex = vertices[begin+1]
        edge_push_intersections!(intersections, prealloc_subgrid,
            first_vertex, second_vertex, nullcline_params, 
            supergrid_dx,
            subgrid_side_n, subgrid_side_x,
            fp_zero_atol
        )
        previous_vertex = second_vertex
        for vertex in vertices[begin+2:end-1]
            interior_push_intersections!(intersections,
                prealloc_subgrid, previous_vertex, vertex, 
                nullcline_params, 
                supergrid_dx,
                subgrid_side_n, subgrid_side_x,
                fp_zero_atol
            )
            previous_vertex = vertex
        end
        edge_push_intersections!(intersections,
            prealloc_subgrid, vertices[end], previous_vertex, # reversed to indicate edge vertex
            nullcline_params, 
            supergrid_dx,
            subgrid_side_n, subgrid_side_x,
            fp_zero_atol
        )
    end
    return collect_array(intersections)
end

_dist(x1::AV, x2::AV) where {AV<:AbstractVector} = sqrt(sum((x1 .- x2) .^ 2)) 
_norm(x1::AV) where {AV<:AbstractVector} = _dist(x1, zero(x1))

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

function interior_push_intersections!(intersections, prealloc_subgrid,
        vertex1, vertex2, nullcline_params, 
        supergrid_dx, subgrid_n, subgrid_x, fp_zero_atol)
    dist_between_vertices = _dist(vertex1, vertex2)
    subgrid_center_dists = 0.:supergrid_dx:dist_between_vertices # step at same res as supergrid
    d_vertex = vertex2 .- vertex1
    d_vertex /= norm(d_vertex)
    for subgrid_center_dist ∈ subgrid_center_dists
        subgrid_center = vertex1 + (subgrid_center_dist * d_vertex)
        interior_subgrid_push_intersections!(intersections,
            prealloc_subgrid, subgrid_center,
            nullcline_params, subgrid_n, subgrid_x,
            fp_zero_atol
        )
    end
end

function calculate_dv_subgrid!(prealloc_subgrid, point::SVector{2,Float64}, params, subgrid_n, subgrid_x)
    subgrid_us = LinRange(max(0., point[1] - subgrid_x/2), min(1., point[1] + subgrid_x/2), subgrid_n) |> collect
    subgrid_vs = LinRange(max(0., point[2] - subgrid_x/2), min(1., point[2] + subgrid_x/2), subgrid_n) |> collect
    @avx for u_idx ∈ axes(subgrid_us, 1), v_idx ∈ axes(subgrid_vs, 1)
        prealloc_subgrid[u_idx, v_idx] = wcm_dv_defn(subgrid_us[u_idx], subgrid_vs[v_idx], params)
    end
    return (subgrid_us, subgrid_vs)
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
  


function interior_subgrid_push_intersections!(intersections, prealloc_subgrid::AbstractArray, point::SVector{2,Float64}, params, subgrid_n::Int, subgrid_x, fp_zero_atol)
    subgrid_us, subgrid_vs = calculate_dv_subgrid!(prealloc_subgrid, point, params, subgrid_n, subgrid_x)
    _interior_subgrid_push_intersections!(intersections, subgrid_us, subgrid_vs, prealloc_subgrid, params, subgrid_x/subgrid_n, fp_zero_atol)
end
function _interior_subgrid_push_intersections!(intersections, subgrid_us, subgrid_vs, subgrid::AbstractArray{<:Number}, params, subgrid_dx, fp_zero_atol)
    # assumes subgrid already calculated
    if !all(sign.(subgrid) .== sign(subgrid[begin]))
        # the sign reverses or is zero on the subgrid, so there is a v_nullcline
        v_nullclines = lines(
            contour(subgrid_us, subgrid_vs, subgrid, 0.)
        )
        _interior_subgrid_push_intersections!(intersections, v_nullclines, params, subgrid_dx, fp_zero_atol)
    end
end

function _interior_subgrid_push_intersections!(intersections, 
        v_nullclines::Vector{<:Curve2}, params, 
        subgrid_dx, fp_zero_atol)
    # if we're here, then a line in v_nullclines must cross a u_nullcline
    for line in v_nullclines
        prev_vertex = SVector{2,Float64}(NaN, NaN)
        prev_du = nothing
        for vertex in line.vertices
            du = wcm_du_defn(vertex..., params)
            if isapprox(du, 0., atol=fp_zero_atol)
                push_intersection!(intersections, vertex, params, subgrid_dx, fp_zero_atol)
            elseif prev_du !== nothing && du * prev_du <= 0
                # crossed zeros
                # find_zero_solution(prev_vertex, vertex, )
                # FIXME save the found solution for future comparison?
                zero_fn = vtx -> wcm_du_defn(vtx..., params)
                interpolated_vertex = find_zero( # just finds du 0 bc in theory vertices define dv=0
                    zero_fn,
                    (prev_vertex, vertex),
                    atol = fp_zero_atol
                )
                # FIXME remove assert
                if abs(wcm_du_defn(interpolated_vertex..., params)) + abs(wcm_dv_defn(interpolated_vertex..., params)) < fp_zero_atol
                    push_intersection!(intersections, interpolated_vertex, params, subgrid_dx, fp_zero_atol)
                end
            end
            prev_vertex = vertex
            prev_du = du
        end
    end
end

function delete_nodes!(node_list)
    for node in node_list
        node.prev.next = node.next
    end
end
du_dv_fn(vtx, params) = [wcm_du_defn(vtx..., params), wcm_dv_defn(vtx..., params)]
function push_intersection!(intersections, intersection_vtx, params, point_dist_atol, fp_zero_atol)
    nearby_existing_intersections = [node for node in intersections if _dist(node.data, intersection_vtx) .< point_dist_atol]
    field_vals = du_dv_fn(intersection_vtx, params)
    sum_abs_field_vals = sum(abs.(field_vals))
    if isempty(nearby_existing_intersections)
        if sum_abs_field_vals < fp_zero_atol
            push!(intersections, intersection_vtx)
        end
    else
        if all(sum(abs.(du_dv_fn(node.data, params))) > sum_abs_field_vals for node in nearby_existing_intersections)
            nearby_existing_intersections[1].data = intersection_vtx
            delete_nodes!(nearby_existing_intersections[begin+1:end])
        end
    end
end



function nt_index(arr::NamedAxisArray{nt_names}, idx) where nt_names
    nt_vals = collect(product(axes_keys(arr)...))[idx]
    NamedTuple{nt_names}(nt_vals)
end
