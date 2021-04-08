calculate_fixedpoints(model::Union{AbstractModel,AbstractSimulation}, args...; kwargs...) = calculate_fixedpoints(get_nullcline_params(model), args...; kwargs...)
calculate_fixedpoints!(du_arr::AbstractArray{T,N}, field_axes, model::Union{AbstractModel{T},AbstractSimulation{T}}, args...; kwargs...) where {N,T} = calculate_fixedpoints!(du_arr, field_axes, get_nullcline_params(model), args...; kwargs...)

function derive_vector_fn(
    nullcline_params::Union{AbstractWCMNullclineParams,AbstractWCMDepNullclineParams}
)
    function wcm_du_dv(uv)
        u = uv[1]; v = uv[2]
        [wcm_du_defn(u, v, nullcline_params)
         wcm_dv_defn(u, v, nullcline_params)]
    end
end

function derive_jacobian_fn(
    nullcline_params::Union{AbstractWCMNullclineParams,AbstractWCMDepNullclineParams}
)
    wcm_du_dv =  derive_vector_fn(nullcline_params)
    uv -> ForwardDiff.jacobian(wcm_du_dv, uv)
end

function derive_vector_fn!(
    nullcline_params::Union{AbstractWCMNullclineParams,AbstractWCMDepNullclineParams}
)
    function wcm_du_dv!(dudv, uv)
        u = uv[1]; v = uv[2]
        dudv[1] = wcm_du_defn(u, v, nullcline_params)
        dudv[2] = wcm_dv_defn(u, v, nullcline_params)
    end
end

function derive_jacobian_fn!(
    nullcline_params::Union{AbstractWCMNullclineParams,AbstractWCMDepNullclineParams}
)
    wcm_du_dv! =  derive_vector_fn!(nullcline_params)
    (result, dudv, uv) -> jacobian!(result, wcm_du_dv!, dudv, uv)
end

function fixedpoint_is_stable(nullcline_params, fp)
    jac_fn = derive_jacobian_fn(nullcline_params)
    jac = jac_fn(fp)
    all(real.(eigvals(jac)) .< 0)
end

function calculate_fixedpoints(
        nullcline_params::Union{AbstractWCMNullclineParams,AbstractWCMDepNullclineParams}, 
        axis_length::Integer=100
        ; 
        kwargs...
    )
    phase_space_bounds = ((0.,1.), (0.,1.))
    field_functions = (wcm_du_defn, wcm_dv_defn)
    calculate_fixedpoints(nullcline_params, field_functions, axis_length, phase_space_bounds; kwargs...)
end

function calculate_fixedpoints(
        nullcline_params::AbstractNullclineParams{T},
        field_functions::NTuple{N,Function},
        axis_length::Integer,
        axes_bounds::NTuple{N,<:Tuple}
        ;
        kwargs...
    ) where {N,T}
    field_axes = [range(b[1], b[2], length=axis_length) for b in axes_bounds]
    field_arr = Array{T,N}(undef, length.(field_axes)...)

    calculate_fixedpoints!(
        field_arr, field_axes, field_functions, 
        nullcline_params, axes_bounds
        ; 
        kwargs...
    )
end

function calculate_fixedpoints!(
        field_arr, field_axes, field_fns, 
        nullcline_params, grand_axes_bounds
        ; 
        step_to_sweep_factor=1.5,
        dist_atol=1e-6, zero_atol=sqrt(eps()),
        kwargs...
    )
    # FIXME? this allocates :(
    nullclines = calculate_nullclines!(
        field_arr, field_axes, field_fns,
        nullcline_params
    )
    potential_intersections = find_potential_intersections!(field_arr, field_axes, nullclines)
    intersection_field_errors = manhattan_norm.(
        apply_field_fns.(Ref(field_fns), potential_intersections, Ref(nullcline_params))
    )
    satisfactory_intersections_idx = intersection_field_errors .< zero_atol
    if all(satisfactory_intersections_idx)
        # could return trivial []
        return potential_intersections
    else
        satisfactory_intersections = FixedPointList( potential_intersections[satisfactory_intersections_idx],
            field_fns, nullcline_params, dist_atol, zero_atol
        ) # FIXME check these hold up to one more zoom in (not hiding another)
        unsatisfactory_intersections = potential_intersections[.!satisfactory_intersections_idx]
        new_satisfactory_intersections_list = recenter_calculate_fixedpoints!.(
            unsatisfactory_intersections,
            Ref(field_arr), Ref(field_axes), Ref(field_fns), 
            Ref(nullcline_params), Ref(grand_axes_bounds);
            step_to_sweep_factor=step_to_sweep_factor,
            dist_atol=dist_atol, zero_atol=zero_atol,
            kwargs...
        )
        for new_satisfactory_intersections ∈ new_satisfactory_intersections_list
            push!.(Ref(satisfactory_intersections), new_satisfactory_intersections)
        end
        return collect_array(satisfactory_intersections)
    end
end

function recenter_calculate_fixedpoints!(new_center,
        field_arr, field_axes, field_fns,
        nullcline_params, grand_axes_bounds; step_to_sweep_factor, kwargs...)
    field_half_width = step(field_axes[begin]) * step_to_sweep_factor / 2
    new_field_axes = [range(max(lb, center_coord - field_half_width), min(ub, center_coord + field_half_width), length=length(field_axes[begin]))
        for ((lb, ub), center_coord) in zip(grand_axes_bounds, new_center)
    ]
    calculate_fixedpoints!(field_arr, new_field_axes, field_fns,
        nullcline_params, grand_axes_bounds
        ;
        step_to_sweep_factor=step_to_sweep_factor, 
        kwargs...
    )
end

### Helper functions (intersections; nullclines; field)

function mark_potential_intersections!(field_arr::AbstractArray{T,N}, field_axes, nullclines) where {T,N}
    axes_length = length(field_axes[begin])
    axes_firsts = first.(field_axes)
    axes_steps = step.(field_axes)
    adjacent_idx_mod = product([[0,1] for _ in 1:N]...)
    contour_count = 0.
    field_arr .= contour_count
    for contour in nullclines
        for line in contour
            for vertex in line.vertices
                # FIXME this will def go out of bounds
                lower_left_vtx = floor.(Int, (vertex .- axes_firsts) ./ axes_steps) .+ 1
                for mod ∈ adjacent_idx_mod
                    idx = lower_left_vtx .+ mod
                    if all(idx .<= axes_length) && field_arr[idx...] == contour_count
                        field_arr[idx...] += 1.
                    end
                end
            end
        end
        contour_count += 1.
    end
    # now potential intersections are marked with contour_count,
    # and everything else has a value < contour_count
end

function find_potential_intersections!(field_arr::AbstractArray{T,N}, field_axes, nullclines) where {T,N}
    if any(isempty.(nullclines))
        return [] # if any nullcline doesn't exist, it can't intersect
    end
    get_ax_vals(idxs::NTuple{N}) = [ax[idx] for (ax, idx) in zip(field_axes, idxs)]
    mark_potential_intersections!(field_arr, field_axes, nullclines)
    intersection_indices = findall(field_arr .== N)
    intersections = [get_ax_vals(Tuple(idxs)) for idxs in intersection_indices]
    return intersections
end

function calculate_nullclines!(
        field_arr, field_axes, field_fns, 
        nullcline_params
    )
    map(field_fns) do field_fn
        calculate_field!(field_arr, field_fn, 
            field_axes, 
            nullcline_params
        )
        lines(contour(field_axes..., field_arr, 0.))
    end
end

function calculate_field!(field_vals::AbstractMatrix{T}, field_fn::Function, usvs::Vector{<:AbstractVector{T}}, p) where T
    us, vs = usvs
    us = collect(us); vs = collect(vs) # FIXME remove line when VectorizationBase updates to v0.18.3(ish)
    @avx for u_idx ∈ axes(us,1), v_idx ∈ axes(vs,1)
        field_vals[u_idx, v_idx] = field_fn(us[u_idx], vs[v_idx], p)
    end
    field_vals
end

function calculate_field(field_fn::Function, usvs::Vector{<:AbstractVector{T}}, p) where T
    us, vs = usvs
    @show us
    dus = Array{T,2}(undef, length(us), length(vs))
    calculate_field!(dus, field_fn, [us,vs], p)
    return dus
end

apply_field_fns(fns, vtx, p) = fns .|> fn -> fn(vtx..., p)

#### Helper structures

manhattan_norm(x) = sum(abs.(x))
euclidean_dist(x1::AV, x2::AV) where {AV<:AbstractVector} = sqrt(sum((x1 .- x2) .^ 2)) 
euclidean_norm(x1::AV) where {AV<:AbstractVector} = _dist(x1, zero(x1))

using DataStructures: MutableLinkedList, ListNode, length
struct FixedPointList{T,MLL<:MutableLinkedList{T},N,FNS<:NTuple{N,Function},P<:AbstractNullclineParams}
    mll::MLL
    field_fns::FNS
    params::P
    dist_atol::Float64
    zero_atol::Float64
    FixedPointList(mll::MLL, field_fns::FNS, params::P, dist_atol, zero_atol) where {T, MLL <: MutableLinkedList{T}, N,FNS<:NTuple{N,Function},P<:AbstractNullclineParams} = new{T,MLL,N,FNS,P}(mll, field_fns, params, dist_atol, zero_atol)
    FixedPointList{T}(field_fns::FNS, params::P, dist_atol, zero_atol) where {T, N,FNS<:NTuple{N,Function},P<:AbstractNullclineParams} = FixedPointList(MutableLinkedList{T}(), field_fns, params, dist_atol, zero_atol)
end
function FixedPointList(points::VEC, args...) where {T,VEC<:AbstractVector{T}}
    fpl = FixedPointList(MutableLinkedList{T}(), args...)
    for point ∈ points
        push!(fpl, point)
    end
    return fpl
end

Base.length(fpl::FixedPointList) = length(fpl.mll)
Base.iterate(l::FixedPointList) = l.mll.len == 0 ? nothing : (l.mll.node.next, l.mll.node.next.next)
Base.iterate(l::FixedPointList, n::ListNode) = n === l.mll.node ? nothing : (n, n.next)
collect_array(l::FixedPointList) = [node.data for node in l]

function delete_nodes!(fps::FixedPointList, node_list::AbstractVector{<:ListNode})
    for node in node_list
        node.prev.next = node.next
        fps.mll.len -= 1
    end
end
function Base.push!(fps::FixedPointList, fixedpoint_vtx)
    field_vals = apply_field_fns(fps.field_fns, fixedpoint_vtx, fps.params)
    new_manhattan_dist = manhattan_norm(field_vals)
    nearby_fps = [node for node in fps if euclidean_dist(node.data, fixedpoint_vtx) < fps.dist_atol]
    if isempty(nearby_fps)
        if new_manhattan_dist < fps.zero_atol
            push!(fps.mll, fixedpoint_vtx)
        end
    else
        if all(manhattan_norm(apply_field_fns(fps.field_fns, node.data, fps.params)) > new_manhattan_dist 
                for node in nearby_fps
            )
            nearby_fps[1].data = fixedpoint_vtx
            delete_nodes!(fps, nearby_fps[begin+1:end])
        end
    end
end