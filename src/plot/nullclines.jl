using WilsonCowanModel: wcm_du_defn, wcm_dv_defn, WCMParams, AbstractNullclineParams, get_nullcline_params
using Makie: @lift

function lifted_wcm_param(;
    Aee=1., Aie=1., Aei=1.5, Aii=0.25,
    θef=0.125, θif=0.4, θib=7., τ=0.4, β=50.,
    decaye=1., decayi=1.)
    @lift WCMParams(;
            Aee=$(Node(Aee)), Aie=$(Node(Aie)), Aei=$(Node(Aei)), Aii=$(Node(Aii)),
            θef=$(Node(θef)), θif=$(Node(θif)), θib=$(Node(θib)), τ=$(Node(τ)), β=$(Node(β)),
            decaye=$(Node(decaye)), decayi=$(Node(decayi))
        )
end

function plot_nullclines!(fig::Figure, p::AbstractNullclineParams, dx=0.01)
    plot_nullclines!(fig, p, 0.:dx:1., 0.:dx:1.)
end
function plot_nullclines!(fig::Figure, p::AbstractNullclineParams, us::AbstractVector, vs::AbstractVector;
        xlabel="u", ylabel="v")
    ax = MakieLayout.Axis(fig, aspect=DataAspect(),
        xlabel=xlabel, ylabel=ylabel
    )

    dus = [wcm_du_defn(u, v, p) for u in us, v in vs]
    dvs = [wcm_dv_defn(u, v, p) for u in us, v in vs]
    u_nullclines = lines(contour(us, vs, dus, 0.))
    v_nullclines = lines(contour(us, vs, dvs, 0.))

    for line in u_nullclines
        xs, ys = coordinates(line)
        Makie.lines!(ax, xs, ys, color=:blue)
    end

    for line in v_nullclines
        xs, ys = coordinates(line)
        Makie.lines!(ax, xs, ys, color=:red, linestyle=:dash, linesize=5)
    end
    return ax
end


function calculate_fixedpoints(mdb_path::AbstractString, new_mods::NamedTuple{NAMES}, dx::T=0.01) where {NAMES, T<:Number}
    prototype_name, saved_mods = read_params_from_data_path(mdb_path)
    @assert NAMES ⊆ keys(saved_mods)
    combined_mods = merge(saved_mods, Dict(pairs(new_mods))) # overwrites saved_mods
    
    variable_mod_names, variable_mod_values = zip([(name, values) for (name,values) in pairs(combined_mods) if length(values)!= 1]...)
    fixedpoints_arr = NamedAxisArray{variable_mod_names}(
        Array{Vector{SVector{2,T}},length(variable_mod_values)}(
            undef, length.(variable_mod_values)...
        ), 
        variable_mod_values...
    )

    prototype = get_prototype(prototype_name)

    map(product(variable_mod_values...)) do mod_values
        fixing_mods = NamedTuple{variable_mod_names}(mod_values)
        # other_opts = Dict() makes sure it saves all frames
        mods = (combined_mods..., fixing_mods..., save_idxs=nothing, save_on=true)
        #write_modifications!(plots_subdir, mods, unique_id)
        model = prototype(; mods...)
        setindex!(fixedpoints_arr, calculate_fixedpoints(model, dx); fixing_mods...)
    end
    return fixedpoints_arr, combined_mods, prototype
end

using Roots
function Roots.find_zero(fn::Function, (p1,p2)::Tuple{PT,PT}, args...) where {N,T, PT<:SVector{N,T}}
    interp(x) = x .* p1 .+ (1-x) .* p2
    interp_zero = find_zero(x -> fn(interp(x)), (0, 1), args...)
    return interp(interp_zero)
end

function calculate_fixedpoints(model::Union{AbstractModel{T},AbstractSimulation{T}}, phase_dx::T=0.01, phase_dy::T=phase_dx) where {T <: Number}
    nullcline_params = get_nullcline_params(model)

    @assert phase_dx == phase_dy
    us = 0.0:phase_dx:1.0
    vs = 0.0:phase_dy:1.0
    radius = sqrt(phase_dx^2 + phase_dy^2)
    
    dus = [wcm_du_defn(u, v, nullcline_params) for u in us, v in vs]
    dvs = [wcm_dv_defn(u, v, nullcline_params) for u in us, v in vs]
    u_nullclines = lines(contour(us, vs, dus, 0.))
    v_nullclines = lines(contour(us, vs, dvs, 0.))

    intersections = []
    prev_point = SVector{2,Float64}(NaN, NaN)
    for u_line in u_nullclines
        prev_dv = 0
        for point in u_line.vertices
            du = wcm_du_defn(point..., nullcline_params)
            dv = wcm_dv_defn(point..., nullcline_params)
            if isapprox(du, 0., atol=eps()) && isapprox(dv, 0., atol=eps())
                push!(intersections, point)
            elseif wcm_dv_defn(point..., nullcline_params) * prev_dv < 0
                # crossed zeros
                # find_zero_solution(prev_point, point, )
                # FIXME save the found solution for future comparison?
                push!(intersections, 
                    find_zero( # just finds dv 0 bc in theory points define du=0
                        pt -> wcm_dv_defn(pt..., nullcline_params), 
                        (prev_point, point)
                    )
                )
            end
            prev_point = point
            prev_dv = dv
        end
    end

    return intersections
end

function nt_index(arr::NamedAxisArray{nt_names}, idx) where nt_names
    nt_vals = collect(product(axes_keys(arr)...))[idx]
    NamedTuple{nt_names}(nt_vals)
end
