
function slice(data::AxisArray{T,2}, target_line::PointVectorLine{2,T,S}, step_fineness=5) where {T,S}
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

function squish(data::AxisArray{T,2}, target_line::PointVectorLine{2,T,S}, args...) where {T,S}
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

# using ForwardDiff
# import ForwardDiff.Dual

# struct DiffCache{T, S}
#     du::Vector{T}
#     dual_du::Vector{S}
# end
# function DiffCache(T, length, ::Type{Val{chunk_size}}) where chunk_size
#     DiffCache(zeros(T, length), zeros(Dual{chunk_size, T}, length))
# end
# DiffCache(u::AbstractArray) = DiffCache(eltype(u),length(u),Val{ForwardDiff.pickchunksize(u |> length)})

function _midpoint(arr::AbstractArray)
    @assert arr[begin] < arr[end]
    return (arr[end] - arr[begin]) / 2 + arr[begin]
end
reduce_along_max_central_gradient(data::NamedAxisArray, args...) = reduce_along_max_central_gradient(data.data, args...)
function reduce_along_max_central_gradient(data::AxisArray{T,2}, 
        reduction::Function=slice, line_fineness=5) where T
    xs, ys = axes_keys(data)
    center_pt = [_midpoint(xs), _midpoint(ys)] 
    interpolation = interpolate(data, Gridded(Linear()))

    neg_gradient_norm(coord) = -norm(Interpolations.gradient(interpolation, coord...))

    result_optim = optimize(neg_gradient_norm, [xs[begin], ys[begin]], 
                                               [xs[end], ys[end]], 
                                               center_pt, 
                                               Fminbox(GradientDescent()); 
                                               )#autodiff=:forward)
    max_grad_coord = Optim.minimizer(result_optim)
    max_grad = Interpolations.gradient(interpolation, max_grad_coord...)

    max_grad_line = try
        PointVectorLine(SA[max_grad_coord...], max_grad)
    catch e
        @warn "Failed to find non-zero max gradient line"
        return (nothing, nothing, nothing, nothing)
    end
    return reduction(interpolation, max_grad_line, xs, ys, line_fineness)
end

reduce_normal_to_halfmax_contour(data::NamedAxisArray, args...) = reduce_normal_to_halfmax_contour(data.data, args...)
function reduce_normal_to_halfmax_contour(data::AxisArray{T,2}, 
        reduction::Function=slice, line_fineness=5) where T
    xs, ys = axes_keys(data)
    center_pt = (_midpoint(xs), _midpoint(ys))
    interpolation = interpolate(data, Gridded(Linear()))

    halfmax_contour_lines = lines(contour(xs, ys, data, 0.5))
    if length(halfmax_contour_lines) == 0
        return (nothing, nothing, nothing, nothing)
    elseif length(halfmax_contour_lines) > 1
       #@warn "Found multiple disconnected contours at level = 0.5: $(length.(coordinates.(halfmax_contour_lines)))"
        return (nothing, nothing, nothing, nothing)

    end
    contour_points = zip(coordinates(only(halfmax_contour_lines))...) |> collect
    nearby = max(div(length(contour_points), 4), 1)
    # if length(contour_points) < 10
    #     @warn "contour very small"
    # end
    cp_center_dx = argmin(map(pt -> norm(pt .- center_pt), contour_points))
    cp_center_pt = contour_points[cp_center_dx]

    nearby_gradients = contour_points[max(1,cp_center_dx-nearby):min(length(contour_points),cp_center_dx+nearby)] .|> pt -> Interpolations.gradient.(Ref(interpolation), pt...)
    resultant_gradient = sum(nearby_gradients)

    resultant_gradient_line = try
        PointVectorLine(SA[cp_center_pt...], resultant_gradient)
    catch e
        @warn "Failed to find non-zero contour-normal gradient"
        return (nothing, nothing, nothing, nothing)
    end
    return reduction(interpolation, resultant_gradient_line, xs, ys, line_fineness)
end

squish(data::NamedAxisArray, args...) = squish(data.data, args...)
slice(data::NamedAxisArray, args...) = slice(data.data, args...)

function draw_reduced_locations!(ax::Nothing, locs)
    @warn "no heatmap provided to inscribe reduction"
end
function draw_reduced_locations!(ax::MakieLayout.Axis, locs)
    xs = [loc[1] for loc in locs]
    ys = [loc[2] for loc in locs]
    plot!(ax, xs, ys, color=:red)
end

function plot_reduction!(scene::Scene, 
                            reducing_fn::Union{typeof(slice), typeof(squish)},
                            data,
                            heatmap_ax::Union{Nothing,MakieLayout.Axis}=nothing)
    reduction_layout = GridLayout()
    dists, vals, locs, line = reduce_normal_to_halfmax_contour(data, reducing_fn)
    if isnothing(dists)
        return nothing
    end
    reduction_ax = MakieLayout.Axis(scene)
    plot!(reduction_ax, dists, vals)
    draw_reduced_locations!(heatmap_ax, locs)
    fitted_sigmoid = fit_sigmoid(vals, dists)
    reduction_title = if fitted_sigmoid != nothing
        sigmoid_val = fitted_sigmoid.(dists)
        plot!(reduction_ax, dists, sigmoid_val, color=:green) 
        Label(scene, 
            "a=$(round(fitted_sigmoid.slope,sigdigits=3)); θ=$(round.(point_from_distance(line, 
                        fitted_sigmoid.threshold),sigdigits=3))", 
            tellwidth=false)
    else
        Label(scene, "no fit", tellwidth=false)
    end
    reduction_ax.xticks = ([dists[begin], dists[end]], 
                                string.([floor.(Ref(Int), locs[begin]), 
                                         floor.(Int, locs[end])]))
    #tightlimits!(reduction_ax)
    ylims!(reduction_ax, 0, 1)
    reduction_layout[:v] = [reduction_title, reduction_ax]
    return reduction_layout
end

# FIXME needs better home
function calc_binary_segmentation(arr)
    never = count(arr .≈ 0)
    always = count(arr .≈ 1)
    total = prod(size(arr))
    sometimes = total - (always + never)
    return (none = never,
            some = sometimes,
            all = always,
            total = total)
end
