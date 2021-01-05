_namedaxisarray_names(naa::NamedAxisArray{names}) where names = names
get_coordinates(aaa::AxisArray) = product(keys.(axes(aaa))...) |> collect


function _collapse_to_axes(A, x_sym::Symbol, y_sym::Symbol)
    name_syms = _namedaxisarray_names(A)
    collapsed_syms = Tuple(setdiff(name_syms, (y_sym, x_sym)))
    data = if findfirst(name_syms .== y_sym) < findfirst(name_syms .== x_sym)
        Simulation73Plotting.avg_across_dims(A, collapsed_syms)'
    else
        Simulation73Plotting.avg_across_dims(A, collapsed_syms)
    end
    return data
end


function _collapse_to_axes(A, x_sym::Symbol)
    name_syms = _namedaxisarray_names(A)
    collapsed_syms = Tuple(setdiff(name_syms, [x_sym]))
    return Simulation73Plotting.avg_across_dims(A, collapsed_syms)
end


get_data(aa::AbstractArray) = aa
get_data(aa::AxisArray) = get_data(parent(aa))

get_axis_values(aa::NamedAxisArray, sym::Symbol) = keys.(axes(aa))[findfirst(AxisIndices.NamedDims.names(aa) .== sym)]
# should be:
# get_axis_values(aa::NamedAxisArray, sym::Symbol) = axes_keys(aa)[dim(aa, sym)]
