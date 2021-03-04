export calculate_fixedpoints!

calculate_fixedpoints(model::Union{AbstractModel,AbstractSimulation}, args...; kwargs...) = calculate_fixedpoints(get_nullcline_params(model), args...; kwargs...)
calculate_fixedpoints!(du_arr::AbstractArray, us::AbstractVector, vs::AbstractVector, model::Union{AbstractModel{T},AbstractSimulation{T}}, args...; kwargs...) where T = calculate_fixedpoints!(du_arr, us, vs, get_nullcline_params(model), args...; kwargs...)

function calculate_fixedpoints(
        nullcline_params::Union{AbstractWCMNullclineParams,AbstractWCMDepNullclineParams}, 
        axis_length::Integer=100
        ; kwargs...
    ) where {T <: Number}
    phase_space_bounds = ((0.,1.), (0.,1.))
    field_functions = (wcm_du_defn, wcm_dv_defn)
    calculate_fixedpoints(nullcline_params, axis_length, phase_space_bounds; kwargs...)
end

function calculate_fixedpoints(
        nullcline_params::AbstractNullclineParams,
        field_functions::NTuple{N,<:Function},
        axis_length::Integer,
        axes_bounds::NTuple{N,<:Tuple};
        kwargs...
    ) where N
    field_axes = [collect(range(b[1], b[2], length=axis_length)) for b in axes_bounds]
    field_arr = Array{T,N}(undef, length.(axes)...)

    calculate_fixedpoints!(
        field_arr, field_axes, field_functions, 
        nullcline_params; 
        kwargs...
    )
end

function calculate_fixedpoints!(
        field_arr, field_axes, field_fns, 
        nullcline_params
        ; 
        manhattan_atol=sqrt(eps())
        kwargs...
    )
    nullclines = calculate_nullclines!(
        field_arr, field_fn, 
        field_axes, 
        nullcline_params
    )
    # if all(.!isempty.(nullclines)), could be intersection
    potential_intersections = find_potential_intersections(FIXME)
    intersection_errors = manhattan_norm.(potential_intersections)
    satisfactory_intersections_idx = intersection_errors .< manhattan_atol
    if all(satisfactory_intersections_idx)
        # could return trivial []
        return potential_intersections
    else
        satisfactory_intersections = potential_intersections[satisfactory_intersections_idx]
        unsatisfactory_intersections = potential_intersections[.!satisfactory_intersections_idx]
        new_satisfactory_intersections = recenter_calculate_fixedpoints!.(
            unsatisfactory_intersections,
            Ref(field_arr), Ref(field_axes), Ref(field_fns), 
            Ref(nullcline_params);
            manhattan_atol = manhattan_atol,
            kwargs...
        )
        push!.(Ref(satisfactory_intersections), new_satisfactory_intersections)
        return satisfactory_intersections
    end
end

manhattan_norm(x) = sum(abs.(x))

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