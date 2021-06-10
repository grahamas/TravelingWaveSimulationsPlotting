using Optim
l2(x,y) = sqrt(sum((x .- y) .^ 2))

function simple_sigmoid(x, a, theta)
    1.0 / (1 + exp(-a * (x - theta)))
end

xs = 0.0:10.0:300.0
ys = simple_sigmoid.(xs, 1.0, 150.0)

slope_est = 1.0
threshold_est = 150.0
lower_param_bounds = [slope_est/5, xs[begin]]
upper_param_bounds = [5*slope_est, xs[end]]

loss(params) = l2(ys, sigmoid(params))

fit_result = optimize(loss, lower_param_bounds, upper_param_bounds,
                          [slope_est, threshold_est], Fminbox(NelderMead()))