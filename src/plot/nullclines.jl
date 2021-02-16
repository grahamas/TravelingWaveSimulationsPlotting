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

# plotting
function nullclines(args...)
    scene, layout = layoutscene()
    layout[1,1] =  nullclines!(scene, args...)
    trim!(layout)
    return (scene, layout)
end

function nullclines!(scene, p::AbstractNullclineParams, dx=0.01)
    nullclines!(scene, p, 0.:dx:1., 0.:dx:1.)
end
function nullclines!(scene, p::AbstractNullclineParams, us::AbstractVector, vs::AbstractVector)

    layout = GridLayout()
    layout[1,1] = ax = MakieLayout.Axis(scene, aspect=DataAspect())

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
    layout[1,0] = Label(scene, "v", tellheight=false)
    layout[end+1,2] = Label(scene, "u", tellwidth=false)
    return layout
end


function calculate_fixedpoints(mdb_path::AbstractString, new_mods::NamedTuple{NAMES}, dx::T=0.01) where {NAMES, T<:Number}
    prototype_name, saved_mods = read_params_from_data_path(mdb_path)
    @assert NAMES ⊆ keys(saved_mods)
    combined_mods = merge(saved_mods, Dict(pairs(new_mods))) # overwrites saved_mods
    
    variable_mod_names, variable_mod_values = zip([(name, values) for (name,values) in pairs(combined_mods) if length(values)!= 1]...)
    n_fixedpoints_arr = NamedAxisArray{variable_mod_names}(zeros(Int, length.(variable_mod_values)...), variable_mod_values...)

    prototype = get_prototype(prototype_name)

    map(product(variable_mod_values...)) do mod_values
        fixing_mods = NamedTuple{variable_mod_names}(mod_values)
        # other_opts = Dict() makes sure it saves all frames
        mods = (combined_mods..., fixing_mods..., save_idxs=nothing, save_on=true)
        #write_modifications!(plots_subdir, mods, unique_id)
        model = prototype(; mods...)
        setindex!(n_fixedpoints_arr, calculate_fixedpoints(model, dx); fixing_mods...)
    end
    return n_fixedpoints_arr, combined_mods, prototype
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

    n_intersections = 0
    prev_point = [NaN, NaN]
    for u_line in u_nullclines
        prev_dv = 0
        for point in u_line.vertices
            du = wcm_du_defn(point..., nullcline_params)
            dv = wcm_dv_defn(point..., nullcline_params)
            if isapprox(du, 0., atol=eps()) && isapprox(dv, 0., atol=eps())
                n_intersections += 1
            elseif wcm_dv_defn(point..., nullcline_params) * prev_dv < 0
                # crossed zeros
                # find_zero_solution(prev_point, point, )
                # FIXME save the found solution for future comparison?
                n_intersections += 1
            end
            prev_point .= point
            prev_dv = dv
        end
    end

    return n_intersections
end

function nt_index(arr::NamedAxisArray{nt_names}, idx) where nt_names
    nt_vals = collect(product(axes_keys(arr)...))[idx]
    NamedTuple{nt_names}(nt_vals)
end
