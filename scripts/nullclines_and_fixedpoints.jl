using TravelingWaveSimulationsPlotting, TravelingWaveSimulations
using WilsonCowanModel, Simulation73
using Dates
using Contour

session_name = "nullclines_and_fixedpoints"
session_id = "$(Dates.now())"

metasweep_A_fpath = TravelingWaveSimulations.get_recent_simulation_data_path(joinpath(homedir(), "data", "ring_normed_blocking", "wider_strength_depthreshold_A"))

function calculate_fixedpoints(mdb_path::AbstractString, new_mods::NamedTuple{NAMES}) where NAMES
    prototype_name, saved_mods = read_params_from_data_path(mdb_path)
    @assert NAMES ⊆ keys(saved_mods)
    combined_mods = merge(saved_mods, Dict(pairs(new_mods))) # overwrites saved_mods
    
    variable_mod_names, variable_mod_values = zip([(name, values) for (name,values) in pairs(combined_mods) if length(values)!= 1]...)
    n_fixedpoints_arr = NamedAxisArray{variable_mod_names}(zeros(Int, length.(variable_mod_values)...), variable_mod_values...)

    map(product(variable_mod_values...)) do mod_values
        fixing_mods = NamedTuple{variable_mod_names}(mod_values)
        # other_opts = Dict() makes sure it saves all frames
        mods = (combined_mods..., fixing_mods..., save_idxs=nothing, save_on=true)
        #write_modifications!(plots_subdir, mods, unique_id)
        prototype = get_prototype(prototype_name)
        model = prototype(; mods...)
        setindex!(n_fixedpoints_arr, calculate_fixedpoints(model); fixing_mods...)
    end
    return n_fixedpoints_arr
end

function calculate_fixedpoints(mdb_path::AbstractString, fixing_mods::Nothing=nothing)
    prototype_name, all_mods = read_params_from_data_path(mdb_path)
    
    fixed_mods = Dict(key => val for (key,val) in pairs(all_mods) if length(val)== 1)
    DUMMY_fixing_mods = Dict(key => val[1] for (key, val) in pairs(all_mods) if length(val) > 1)

    # other_opts = Dict() makes sure it saves all frames
    mods = (fixed_mods..., DUMMY_fixing_mods..., save_idxs=nothing, save_on=true)

    #write_modifications!(plots_subdir, mods, unique_id)
    prototype = get_prototype(prototype_name)
    model = prototype(; mods...)

    return translate_model_to_nullcline_params(model)

    return calculate_fixedpoints(model)
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

function calculate_fixedpoints(model::Union{AbstractModel{T},AbstractSimulation{T}}) where T
    nullcline_params = translate_model_to_nullcline_params(model)
    
    dx = dy = 0.1
    us = 0.0:dx:1.0
    vs = 0.0:dy:1.0
    radius = sqrt(dx^2 + dy^2)
    
    dus = [WilsonCowanModel.wcm_du_defn(u, v, nullcline_params) for u in us, v in vs]
    dvs = [WilsonCowanModel.wcm_dv_defn(u, v, nullcline_params) for u in us, v in vs]
    u_nullclines = Contour.lines(Contour.contour(us, vs, dus, 0.))
    v_nullclines = Contour.lines(Contour.contour(us, vs, dvs, 0.))

    n_intersections = 0
    for u_line in u_nullclines
        for v_line in v_nullclines
            n_intersections += count_intersections(u_line, v_line, dx)
        end
    end

    return n_intersections
end

calc_distance((p1, p2)) = sqrt(sum((p1 .- p2) .^ 2))
function count_intersections(line1, line2, radius)
    # 1. find local minima of distance between vertices
    # 2. stupid thing: if distance <= radius (sqrt(dx^2 + dy^2)) then they cross
    # FIXME account for parallel lines?
    # FIXME there's a faster way to do this

    point_pairs = product(line1.vertices, line2.vertices)
    intersections = findall(calc_distance.(point_pairs) .<= radius)
    return length(intersections)
end

calculate_fixedpoints(metasweep_A_fpath, (stim_strength=2.0, blocking_θI = 9.,
    #Aii=125., Aie=125.)
    Aii=1.0, Aie=15.0)
)