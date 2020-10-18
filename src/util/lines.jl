l2(x,y) = sqrt(sum((x .- y) .^ 2))
function diagonal_slice(x_axis, y_axis, data::Matrix, y_intercept, slope, dx=1)
    interp = LinearInterpolation((x_axis, y_axis), data')
    calc_y(x) = slope * x + y_intercept
    sample_points = [(x, calc_y(x)) for x in x_axis if y_axis[begin] <= calc_y(x) <= y_axis[end]]
    sample_values = sample_points .|> (x) -> interp(x...)
    distance_along_line = l2.(Ref(sample_points[begin]), sample_points)
    return (distance_along_line, sample_points, sample_values)
end

function project(coord::C, line::C, line_point::C) where {N, C <: SVector{N}}
    projection = (line * line') / (line' * line)
    return projection * coord + (I - projection) * line_point
end

line_dist_to_coord(dist, line, origin) = (dist * line) + origin 
#axes_vals(data::AbstractAxisArray) = keys.(axes(data))

struct PointVectorLine{N,T,S<:SVector{N,T}}
    point::S
    vector::S  # FIXME should have special vertical and horizontal lines and/or insist norm(vector) != 0
    PointVectorLine(point::S,vector::S) where {N,T,S<:SVector{N,T}} = begin
        @assert any(abs.(vector) .> sqrt(eps(T)))
        new{N,T,S}(point,vector)
    end
end
function slope(line::PointVectorLine{2})
    if any(line.vector .== 0.0) # FIXME can't be vertical; shouldn't be allowed
        return 0.0
    else
        return line.vector[2] / line.vector[1]
    end
end
# the point where line takes value val in dim dimension
function point_from_dim_val(line::PointVectorLine, val::Number, dim::Int)
    scale = (val - line.point[dim]) / line.vector[dim]
    return line.point .+ line.vector .* scale
end
x_from_y(line::PointVectorLine{2}, y::Number) = point_from_dim_val(line, y, 2)[1]
y_from_x(line::PointVectorLine{2}, x::Number) = point_from_dim_val(line, x, 1)[2]
y_intercept(line::PointVectorLine{2}) = y_from_x(line, 0.0)
x_intercept(line::PointVectorLine{2}) = x_from_y(line, 0.0) 
point_from_distance(line::PointVectorLine, dist::Number) = line.vector .* dist .+ line.point
point_from_distance(::Any, ::Missing) = missing
get_orthogonal_vector(line::PointVectorLine) = -SA[line.vector[2], -line.vector[1]]

function originate_from_left(line, xs, ys)
    
    x_min, x_max = extrema(xs)
    y_min, y_max = extrema(ys)
    if all(line.vector .== 0.0)
        return PointVectorLine(SA[x_min, line.point[2]], SA[1.0, 0.0])
    elseif line.vector[1] == 0.0
        return PointVectorLine(SA[line.point[1], y_min], SA[0.0, 1.0])
    end

    # need increasing x
    new_vector = line.vector[1] > 0 ? line.vector : -line.vector

    if y_min <= y_from_x(line, x_min) <= y_max
        return PointVectorLine(SA[x_min, y_from_x(line, x_min)], new_vector)
    elseif x_min <= x_from_y(line, y_max) <= x_from_y(line, y_min)
        @assert x_from_y(line, y_max) <= x_max
        return PointVectorLine(SA[x_from_y(line, y_max), y_max], new_vector)
    elseif x_min <= x_from_y(line, y_min) <= x_from_y(line, y_max)
        @assert x_from_y(line, y_min) <= x_max
        return PointVectorLine(SA[x_from_y(line, y_min), y_min], new_vector)
    else
        @show line
        @show x_from_y(line, y_min)
        @show x_from_y(line, y_max)
        @show y_from_x(line, x_min)
        @show y_from_x(line, x_max)
        @show extrema(xs)
        @show extrema(ys)
        error("Line does not fall within axes")
    end
end