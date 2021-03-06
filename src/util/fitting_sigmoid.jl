struct FittedSigmoid{T}
    left_val::T
    change::T
    slope::T
    threshold::T
    error::T
end
FittedSigmoid(::Missing) = FittedSigmoid{Missing}(missing, missing, missing, missing, missing)

(s::FittedSigmoid)(x) = s.left_val .+ s.change .* NeuralModels.simple_sigmoid(x, s.slope, s.threshold)

fit_sigmoid(::Nothing, ::Nothing) = FittedSigmoid(missing)
function fit_sigmoid(ys, xs)
    # fits sigmoid that goes from 0 to 1
    if xs[end] < xs[begin]
        xs = reverse(xs)
        ys = reverse(ys)
    end
    if xs[begin] == xs[end] || length(xs) < 3
        return FittedSigmoid(missing)
    end
    left_val = ys[begin]
    sigmoid_change = ys[end] - ys[begin]
    sigmoid(params) = left_val .+ sigmoid_change .* NeuralModels.simple_sigmoid.(xs, params[1], params[2])
    loss(params) = l2(ys, sigmoid(params))
    D = CenteredDifference{1}(1, 2, xs[2] - xs[1], length(ys))
    derivative_estimate = D*ys
    slope_est, i_max = findmax(abs.(derivative_estimate[begin+1:end-1]))
    i_max += 1 # ensure i_max is not on the boundary
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

# offset - A sigmoid(-a (x - theta))
# where x is distance along PointVectorLine(θ, φ)
struct SigmoidOnSlice{T}
    offset::T
    A::T
    a::T
    θ::SVector{2,T}
    φ::SVector{2,T}
    error::T
end
function SigmoidOnSlice(sig::FittedSigmoid{T}, lin::PointVectorLine) where {T<:Number}
    SigmoidOnSlice{T}(
        sig.left_val,
        sig.change,
        sig.slope,
        point_from_distance(lin, sig.threshold),
        lin.vector,
        sig.error
    )
end
function SigmoidOnSlice(sig::FittedSigmoid{Missing}, ::Any)
    SigmoidOnSlice{Missing}(
        missing,
        missing,
        missing,
        SA[missing, missing],
        SA[missing, missing],
        missing        
    )
end