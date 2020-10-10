using TravelingWaveSimulationsPlotting: reduce_normal_to_halfmax_contour

using IterTools, AxisIndices, Interpolations, StaticArrays

function simple_sigmoid_fn(x, a, theta)
    1.0 / (1 + exp(-a * (x - theta)))
end

xs, ys = 0.0:10.0:300.0, 0.0:10.0:290.0
grid = product(xs, ys)

sigmoid_unit_diagonal(x, y) = simple_sigmoid_fn(x - y, 0.1, 0.0)

field_unit_diagonal = AxisArray([1 - sigmoid_unit_diagonal(x,y) + 0.1 * rand() for (x,y) in grid], xs, ys);

using Makie: heatmap
plot_sigmoid_unit_diagonal = heatmap(xs, ys, field_unit_diagonal.parent);

dists, vals, locs, line = reduce_normal_to_halfmax_contour(field_unit_diagonal)

using Makie: scatter!
scatter!(plot_sigmoid_unit_diagonal, locs)

