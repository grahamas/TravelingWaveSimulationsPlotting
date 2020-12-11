_namedaxisarray_names(naa::NamedAxisArray{names}) where names = names
get_coordinates(aaa::AxisArray) = product(keys.(axes(aaa))...) |> collect


function _collapse_to_axes(A, x_sym, y_sym)
    name_syms = _namedaxisarray_names(A)
    collapsed_syms = Tuple(setdiff(name_syms, (y_sym, x_sym)))
    data = if findfirst(name_syms .== y_sym) < findfirst(name_syms .== x_sym)
        Simulation73Plotting.avg_across_dims(A, collapsed_syms)'
    else
        Simulation73Plotting.avg_across_dims(A, collapsed_syms)
    end
    return data
end

get_data(aa::AbstractArray) = aa
get_data(aa::AxisArray) = get_data(parent(aa))
