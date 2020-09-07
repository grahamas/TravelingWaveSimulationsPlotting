
function slice(data::AbstractAxisArray{T,2}, target_line::PointVectorLine{2,T,S}, step_fineness=5) where {T,S}
    xs, ys = axes_keys(data)
    interpolation = interpolate(data, Gridded(Linear()))
    return slice(interpolation, target_line, xs, ys, step_fineness)
end

function slice(interpolation::AbstractInterpolation, target_line::PointVectorLine{2,T,S}, xs, ys, step_fineness=5) where {T,S}
    # We want this line to proceed from its point through the grid (left to right for convention's and plotting's sake) 
    line = originate_from_left(target_line, xs, ys)

    # Calculate a reasonable step
    dx = abs(xs[begin+1] - xs[begin])
    dy = abs(xs[begin+1] - xs[begin])
    step = sqrt(dx^2 + dy^2) / step_fineness

    # Step along line until
    distance_along_line = 0
    dists = T[]; points = S[]; values = T[];
    next_point = point_from_distance(line, distance_along_line)
    while xs[begin] <= next_point[1] <= xs[end] && ys[begin] <= next_point[2] <= ys[end]
        next_value = interpolation(next_point...)
        # FIXME preallocate
        push!(dists, distance_along_line)
        push!(points, next_point)
        push!(values, next_value)
        distance_along_line += step
        next_point = point_from_distance(line, distance_along_line)
    end
    return (dists, values, points, line) 
end

function squish(data::AbstractAxisArray{T,2}, target_line::PointVectorLine{2,T,S}, args...) where {T,S}
    xs, ys = axes_keys(data)
    interpolation = interpolate(data, Gridded(Linear()))
    return squish(interpolation, target_line, xs, ys, args...)
end
function squish(interpolation::AbstractInterpolation, target_line::PointVectorLine{2,T,S}, xs, ys, line_fineness=5, squish_fineness=2) where {T,S}
    squished_dists, _, squished_points, squished_line = slice(interpolation, target_line, xs, ys, line_fineness)
    orthogonal_vec = get_orthogonal_vector(target_line)
    squished_vals = map(squished_points) do point
        _, crosscut_vals, _ = slice(interpolation, PointVectorLine(point, orthogonal_vec), xs, ys, squish_fineness)
        mean(crosscut_vals)
    end
    return squished_dists, squished_vals, squished_points, squished_line
end

function _midpoint(arr::AbstractArray)
    @assert arr[begin] < arr[end]
    return (arr[end] - arr[begin]) / 2 + arr[begin]
end
function reduce_along_max_central_gradient(data::AbstractAxisArray{T,2}, reduction::Function=slice, line_fineness=5) where T
    xs, ys = axes_keys(data)
    center_pt = [_midpoint(xs), _midpoint(ys)] 
    interpolation = interpolate(data, Gridded(Linear()))
    neg_gradient_norm(coord) = -norm(Interpolations.gradient(interpolation, coord...))
    result_optim = optimize(neg_gradient_norm, [xs[begin], ys[begin]], 
                                               [xs[end], ys[end]], 
                                               center_pt, 
                                               Fminbox(GradientDescent()); 
                                               autodiff=:forward)
    max_grad_coord = Optim.minimizer(result_optim)
    max_grad = Interpolations.gradient(interpolation, max_grad_coord...)
    max_grad_line = PointVectorLine(SA[max_grad_coord...], max_grad)
    @show string(reduction)
    return reduction(interpolation, max_grad_line, xs, ys, line_fineness)
end

squish(data::NamedAxisArray, args...) = squish(data.data, args...)
slice(data::NamedAxisArray, args...) = slice(data.data, args...)
reduce_along_max_central_gradient(data::NamedAxisArray, args...) = reduce_along_max_central_gradient(data.data, args...)

function draw_reduced_locations!(ax::Nothing, locs)
    @warn "no heatmap provided to inscribe reduction"
end
function draw_reduced_locations!(ax::LAxis, locs)
    xs = [loc[1] for loc in locs]
    ys = [loc[2] for loc in locs]
    plot!(ax, xs, ys, color=:red)
end

function plot_max_gradient!(scene::Scene, 
                            reducing_fn::Union{typeof(slice), typeof(squish)},
                            data,
                            heatmap_ax::Union{Nothing,LAxis}=nothing)
    max_grad_layout = GridLayout()
    dists, vals, locs, line = reduce_along_max_central_gradient(data, slice)
    max_grad_ax = LAxis(scene)
    plot!(max_grad_ax, dists, vals)
    draw_reduced_locations!(heatmap_ax, locs)
    fitted_sigmoid = fit_sigmoid(vals, dists)
    max_grad_title = if fitted_sigmoid != nothing
        sigmoid_val = fitted_sigmoid.(dists)
        plot!(max_grad_ax, dists, sigmoid_val, color=:green) 
        LText(scene, 
            "a=$(round(fitted_sigmoid.slope,sigdigits=3)); θ=$(round.(point_from_distance(line, 
                        fitted_sigmoid.threshold),sigdigits=3))", 
            tellwidth=false)
    else
        LText(scene, "no fit", tellwidth=false)
    end
    max_grad_ax.xticks = ([dists[begin], dists[end]], 
                                string.([floor.(Ref(Int), locs[begin]), 
                                         floor.(Int, locs[end])]))
    #tightlimits!(max_grad_ax)
    ylims!(max_grad_ax, 0, 1)
    max_grad_layout[:v] = [max_grad_title, max_grad_ax]
    return max_grad_layout
end

# FIXME needs better home
function calc_binary_segmentation(arr)
    never = sum(arr .== 0)
    always = sum(arr .== 1)
    total = prod(size(arr))
    sometimes = total - (always + never)
    return (none = never / total,
            some = sometimes / total,
            all = always / total)
end
