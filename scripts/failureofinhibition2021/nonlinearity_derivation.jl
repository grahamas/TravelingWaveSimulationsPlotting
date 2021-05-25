#using GLMakie; ext_2d = "png"; GLMakie.activate!()
using CairoMakie; ext_2d = "svg"; CairoMakie.activate!(); 
using DrWatson, Dates
using Colors

let figure_name = "nonlinearity_derivation",
    N = 1000,
    fig_resolution=(800,600),
    simple_theme = Theme(
        linewidth = 3.0,
        fontsize=24,
        backgroundcolor =  RGBA(1.,1.,1.,0.),
        Axis = (
            backgroundcolor = RGBA(1.,1.,1.,0.),
            leftspinevisible = true,
            rightspinevisible = false,
            bottomspinevisible = true,
            topspinevisible = false,
            xgridcolor = RGBA(1.,1.,1.,0.),
            ygridcolor = RGBA(1.,1.,1.,0.),
        )
    ) ;

_subplotsdir = plotsdir("$(figure_name)_$(ext_2d)_$(now())")
subplotsdir(x="") = joinpath(_subplotsdir, x)

with_theme(simple_theme) do
θ_firing = randn(N) .* 0.01 .+ 0.1
θ_failing = randn(N) .* 0.01 .+ 0.2

is_firing(x)= x ? 1. : 0.

example_firing_threshold = 0.1
example_failing_threshold = 0.2

xs = collect(0.0:0.001:0.3)

example_monotonic_nonl = is_firing.(example_firing_threshold .<= xs)
example_failing_nonl = is_firing.(example_firing_threshold .<= xs .<= example_failing_threshold)

fig = Figure(resolution=fig_resolution)
fig[1,1] = ax_firing_example = Axis(fig, ylabel="firing?",
    yticks=[0,1])
lines!(ax_firing_example, xs, example_monotonic_nonl)

firing_accum = zeros(size(xs))
fig[2,1] = ax_firing_many = Axis(fig, ylabel="firing?",
    yticks=[0,1])
for θ ∈ θ_firing
    this_firing = is_firing.(θ .<= xs)
    lines!(ax_firing_many, xs, this_firing)
    firing_accum .+= this_firing
end
firing_accum ./= length(θ_firing)
fig[3,1] = ax_firing_accum = Axis(fig,
    xlabel="input", ylabel="prop. firing")
lines!(ax_firing_accum, xs, firing_accum)

fig[1,2] = ax_failing_example = Axis(fig)
lines!(ax_failing_example, xs, example_failing_nonl)

failing_accum = zeros(size(xs))
fig[2,2] = ax_failing_many = Axis(fig)
for (lo, hi) ∈ zip(θ_firing, θ_failing)
    this_failing = is_firing.(lo .<= xs .<= hi)
    lines!(ax_failing_many, xs, this_failing)
    failing_accum .+= this_failing
end
failing_accum ./= length(θ_failing)
fig[3,2] = ax_failing_accum = Axis(fig)
lines!(ax_failing_accum, xs, failing_accum)

xlims!.([ax_firing_example, ax_failing_example,
         ax_firing_many, ax_failing_many,
         ax_firing_accum, ax_failing_accum], Ref((xs[begin], xs[end])))
ylims!.([ax_firing_example, ax_failing_example,
         ax_firing_many, ax_failing_many,
         ax_firing_accum, ax_failing_accum], Ref((0., 1.01)))
hidedecorations!.([ax_failing_example,
                   ax_failing_many,
                   ax_failing_accum])
hidexdecorations!.([ax_firing_many, ax_firing_example])

fig[0,1] = Label(fig, "fire", tellwidth=false)
fig[1,2] = Label(fig, "fire → fail", tellwidth=false)

mkpath(subplotsdir())
save(subplotsdir("$(figure_name).$(ext_2d)"), fig)

fig

end # with_theme
end # let