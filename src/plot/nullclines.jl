using WilsonCowanModel: wcm_du_defn, wcm_dv_defn, WCMParams, AbstractNullclineParams, wcm_nullcline_params
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
    return (scene, layout)
end

function nullclines!(scene, p::AbstractNullclineParams, dx=0.01)
    nullclines!(scene, p, 0.:dx:1., 0.:dx:1.)
end
function nullclines!(scene, p::AbstractNullclineParams, us::AbstractVector, vs::AbstractVector)

    layout = GridLayout()
    layout[1,1] = ax = LAxis(scene)

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
    layout[1,0] = LText(scene, "v", tellheight=false)
    layout[end+1,2] = LText(scene, "u", tellwidth=false)
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

using ForwardDiff

function objective_away_from_roots(fn, cache, roots=[])
    f(u) = begin 
        F = WilsonCowanModel.wcm(u, nullcline_params, 0.)
        root_dists = sum.((F .- roots) .^ 2)
        root_closenesses = 1 ./ root_dists
        total_root_closeness = sum(root_closenesses)
        sum(F .^ 2) * total_root_closeness
    end
    g!(result, u) = ForwardDiff.gradient!(result, f, u)
    h!(result, u) = ForwardDiff.hessian!(result, f, u)
    return TwiceDifferentiable(f, g!, h!, cache)
end

# function calculate_fixedpoints(model::Union{AbstractModel{T},AbstractSimulation{T}}) where T
#     nullcline_params = translate_model_to_nullcline_params(model)
#     # f(u::Vector{T}) where T = WilsonCowanModel.wcm(u, nullcline_params, 0.)
#     # g(u::Vector{T}) where T = ForwardDiff.jacobian(f, u)
#     find_all_roots(objective, ([0., 0.], [1., 1.]))
# end

using Optim, IterTools
function find_all_roots(fn::Function, (lower, upper)::Tuple{Vector{T}, Vector{T}}, roots=[], max_recursion=8) where T
    if max_recursion <= 0
        @warn "find_all_roots recursed too far"
        return []
    end
    objective = objective_away_from_roots(
        u -> WilsonCowanModel.wcm(u, nullcline_params, 0.), 
        zeros(T, 2), 
        roots
    )
    x0 = (upper .- lower) ./ 2. .+ lower
    @show lower x0 upper
    constraints = TwiceDifferentiableConstraints(lower, upper)
    result = optimize(objective, constraints, x0, IPNewton())
    @show Optim.minimum(result)
    if Optim.converged(result) && isapprox(Optim.minimum(result), 0., atol=sqrt(eps()))
        root = Optim.minimizer(result)
        lower_coords = product(zip(lower, root .+ sqrt(eps(T)),)...) .|> collect
        upper_coords = product(zip(root .- sqrt(eps(T)), upper)...) .|> collect
        quadrant_bounds = zip(lower_coords, upper_coords) |> collect
        return [root for root in find_all_roots.(f, g!, h!, quadrant_bounds, max_recursion-1) if !isempty(root)] # FIXME could avoid allocations if we pass objective around
    else
        return []
    end
end

struct Intersection{PointType,T}
    identified_points::Vector{PointType}
    error::T
end
struct Intersections{PointType,T}
    intersections::Vector{Intersection{PointType,T}}
    error::T
end
function add_point!(intersection::Intersection{P}, point::P) where P
    if any(calc_distance.(intersection.identified_points, Ref(point)) .<= intersection.error)
        push!(intersection.identified_points, point)
        return true
    else 
        return false
    end
end

function add_putative_intersection!(radius::Number, point::P) where P
    Intersections([Intersection([point], radius)], radius)
end
function add_putative_intersection!(existing_intersections::Intersections{P}, point::P) where P
    """
        Adds point to list of intersections, either:
            1) as new intersection
            2) as "same" as existing intersection
            3) as "same" as multiple existing intersections, combining
    """
    intersections_arr = existing_intersections.intersections
    impinged_intersections = add_point!.(intersections_arr, Ref(point))
    if sum(impinged_intersections) == 0
        # make new intersection
        push!(intersections_arr, Intersection([point], existing_intersections.error))
    elseif sum(impinged_intersections) > 1
        # some intersections can be joined
        keep_idx, dispose_idxs... = findall(impinged_intersections)
        for dispose_idx in dispose_idxs
            append!(intersections_arr[keep_idx].identified_points, intersections_arr[dispose_idx].identified_points)
        end
        deleteat!(intersections_arr, dispose_idxs)
    end 
    return existing_intersections
end

# function add_putative_intersection!(existing_intersections::Intersections{P}, point::P) where P
#     push!(existing_intersections.intersections, Intersection([point], existing_intersections.error))
#     return existing_intersections
# end


function calculate_fixedpoints(model::Union{AbstractModel{T},AbstractSimulation{T}}, phase_dx::T=0.01, phase_dy::T=phase_dx) where {T <: Number}
    nullcline_params = wcm_nullcline_params(model)

    @assert phase_dx == phase_dy
    us = 0.0:phase_dx:1.0
    vs = 0.0:phase_dy:1.0
    radius = sqrt(phase_dx^2 + phase_dy^2)
    
    dus = [wcm_du_defn(u, v, nullcline_params) for u in us, v in vs]
    dvs = [wcm_dv_defn(u, v, nullcline_params) for u in us, v in vs]
    u_nullclines = lines(contour(us, vs, dus, 0.))
    v_nullclines = lines(contour(us, vs, dvs, 0.))

    n_intersections = 0
    for u_line in u_nullclines
        for v_line in v_nullclines
            n_intersections += count_intersections(u_line, v_line, radius)
        end
    end

    return n_intersections
end

calc_distance(p1::P, p2::P) where P = calc_distance((p1,p2))
calc_distance((p1, p2)) = sqrt(sum((p1 .- p2) .^ 2))
function count_intersections(line1, line2, radius)
    # 1. find local minima of distance between vertices
    # 2. stupid thing: if distance <= radius (sqrt(dx^2 + dy^2)) then they cross
    # FIXME account for parallel lines?
    # FIXME there's a faster way to do this

    all_point_pairs = product(line1.vertices, line2.vertices)
    adjacent_pair_idxs = findall(calc_distance.(all_point_pairs) .<= radius)
    length(adjacent_pair_idxs) == 0 && return 0
    intersections = radius
    for adjacent_pair_idx in adjacent_pair_idxs
        intersections = check_adjacent_pair_as_intersection!(intersections,
            line1.vertices, line2.vertices, adjacent_pair_idx)
    end
    intersections isa Number && return 0
    return length(intersections.intersections)
end



function nt_index(arr::NamedAxisArray{nt_names}, idx) where nt_names
    nt_vals = collect(product(axes_keys(arr)...))[idx]
    NamedTuple{nt_names}(nt_vals)
end

function check_adjacent_pair_as_intersection!(intersections, vertices1, vertices2, pair_idx)
    v1_idx, v2_idx = Tuple(pair_idx)
    seg1_idxs = max(1,v1_idx-1):min(length(vertices1),v1_idx+1)
    seg2_idxs = max(1,v2_idx-1):min(length(vertices2),v2_idx+1)
    seg1 = vertices1[seg1_idxs]
    seg2 = vertices2[seg2_idxs]
    if segments_cross(seg1, seg2)
        return add_putative_intersection!(intersections, vertices1[v1_idx])
    else
        return intersections
    end
end

function segments_cross(seg1, seg2, atol = eps())
    # first check if any endpoints are within atol
    any(map(zip(seg1, seg2)) do (p1, p2)
        # true if points equal
        return all(isapprox.(p1, p2, atol=atol))
    end) && return true  # two endpoints approx equal

    s1p1, s1p2 = seg1
    s2p1, s2p2 = seg2
    p1_polarities = s1p1 .- s2p1
    p2_polarities = s1p2 .- s2p2

    # any polarity reverses <=> segments crossed
    return any(p1_polarities .* p2_polarities .< 0)
end

