struct FittedSigmoid{T}
    left_val::T
    change::T
    slope::T
    threshold::T
    error::T
end

(s::FittedSigmoid)(x) = s.left_val .+ s.change .* NeuralModels.simple_sigmoid_fn(x, s.slope, s.threshold)

function fit_sigmoid(ys, xs)
    # fits sigmoid that goes from 0 to 1
    if xs[end] < xs[begin]
        xs = reverse(xs)
        ys = reverse(ys)
    end
    if xs[begin] == xs[end]
        return nothing
    end
    left_val = ys[begin]
    sigmoid_change = ys[end] - ys[begin]
    sigmoid_fn(params) = left_val .+ sigmoid_change .* NeuralModels.simple_sigmoid_fn.(xs, params[1], params[2])
    loss(params) = l2(ys, sigmoid_fn(params))
    D = CenteredDifference{1}(1, 2, xs[2] - xs[1], length(ys))
    derivative_estimate = D*ys
    slope_est, i_max = findmax(abs.(derivative_estimate))
    threshold_est = xs[i_max]

    lower_param_bounds = [slope_est/5, xs[begin]]
    upper_param_bounds = [5*slope_est, xs[end]]
    fit_result = optimize(loss, lower_param_bounds, upper_param_bounds,
                          [slope_est, threshold_est], Fminbox(NelderMead());
                          autodiff=:forward)
    slope_fit, threshold_fit = Optim.minimizer(fit_result)
    error_fit = Optim.minimum(fit_result) |> abs
    return FittedSigmoid(
                         left_val,
                         sigmoid_change,
                         slope_fit,
                         threshold_fit,
                         error_fit
                        )
end

function collapse_and_fit_sigmoid(uncollapsed_space, (x_axis_name, y_axis_name), reduction::Function=squish)
    phase_space = _collapse_to_axes(uncollapsed_space, x_axis_name, y_axis_name)
    line_dists, line_vals, line_locs, line = reduce_along_max_central_gradient(phase_space, reduction)
    fit_sigmoid(line_vals, line_dists)
end