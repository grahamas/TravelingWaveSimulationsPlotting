function _collapse_to_axes(A::NamedDimsArray{Names}, x_sym::Symbol, y_sym::Symbol) where Names
    collapsed_syms = Tuple(setdiff(Names, (y_sym, x_sym)))
    data = if findfirst(Names .== y_sym) < findfirst(Names .== x_sym)
        Simulation73Plotting.avg_across_dims(A, collapsed_syms)'
    else
        Simulation73Plotting.avg_across_dims(A, collapsed_syms)
    end
    return data
end


function _collapse_to_axes(A::NamedDimsArray{Names}, x_sym::Symbol) where Names
    collapsed_syms = Tuple(setdiff(Names, [x_sym]))
    return Simulation73Plotting.avg_across_dims(A, collapsed_syms)
end
