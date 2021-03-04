using Dates#, LaTeXStrings

const POPS = [:E,:I]
LABEL_DICT = merge(
    Dict(Symbol(var,pop1,pop2) => "$(var)_{$(pop1)$(pop2)}" for 
        var=[:A,:S], pop1=POPS, pop2=POPS
    )
)

function default_label_translate(sym::Symbol)
    getkey(LABEL_DICT, sym, sym |> string)
end

NaN_if_missing(::Val{T}, ::Missing) where {T} = T(NaN)
NaN_if_missing(::Val{T}, val::T) where {T} = val


function _save!(fig, plot_name;
        unique_id="$(Dates.now())",
        session_id=unique_id,
        session_name="unnnamedsession",
        plots_subdir=plotsdir("$(session_name)_$(session_id)"))
    mkpath(plots_subdir)
    save_path = plotsdir(plots_subdir, plot_name)
    @info "saving $save_path"
    Makie.save(save_path, fig)
    return fig
end


function plot_and_save(plot_fn!, args...; 
        unique_id="$(Dates.now())",
        session_id=unique_id,
        session_name="unnnamedsession",
        plots_subdir=plotsdir("$(session_name)_$(session_id)"),
        plot_name = "$(strip(string(plot_fn!), '!'))_$(unique_id).png",
        figure_resolution=(1600, 1600),
        kwargs...)
    figure = Figure(resolution=figure_resolution)

    figure[1,1] = plot_fn!(figure, args...; kwargs...)

    _save!(figure, plot_name; unique_id=unique_id, session_id=session_id, session_name=session_name, plots_subdir=plots_subdir)

    return figure
end

export plot_and_save_many
function plot_and_save_many(plot_fn, args...; 
        unique_id="$(Dates.now())",
        plot_name="$(strip(string(plot_fn), '!'))_$(unique_id).png",
        session_name="unnnamedsession",
        session_id=unique_id,
        kwargs...)
    "Like plot_and_save, but assuming plot_fn! returns many figures"

    figs = map(plot_fn(args...; kwargs...)) do (fig, layout, name)
        _save!(fig, "$(name)_$(plot_name)";
            session_id=session_id,
            session_name=session_name)
        fig
    end
    return figs
end

function layout_plot(plot_fn, args...; kwargs...)
    fig = Figure()
    fig[1,1] =  plot_fn(fig, args...; kwargs...)
    return (fig, layout)
end

function figure_plot(plot_fn, args...; resolution = nothing, kwargs...)
    figure = resolution === nothing ? Figure() : Figure(resolution = resolution)
    figure[1,1] =  plot_fn(figure, args...; kwargs...)
    return figure
end