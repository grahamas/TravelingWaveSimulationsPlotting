using NeuralModels: simple_sigmoid_fn
using GLMakie

my_model(x, (a_1, a_2), (θ_1, θ_2)) = simple_sigmoid_fn(x, a_1, θ_1) - simple_sigmoid_fn(x, a_2, θ_2)

kim_model(x, (a_1, a_2), (θ_1, θ_2)) = simple_sigmoid_fn(x, a_1, θ_1)  * (1-simple_sigmoid_fn(x, a_2, θ_2))

low_gain = (1.0, 1.0)
high_gain = (50.0, 50.0)

close_thresholds = (1.0, 2.0)
far_thresholds = (1.0, 7.0)

xs = -3.:0.01:10.

fig = Figure()
fig[1,1] = ax = Axis(fig)
lines!(ax, xs, my_model.(xs, Ref(low_gain), Ref(close_thresholds)))
lines!(ax, xs, kim_model.(xs, Ref(low_gain), Ref(close_thresholds)))