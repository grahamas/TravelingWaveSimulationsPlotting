export fixedpoint_is_stable, fixedpoint_stability, 
    fixedpoint_is_oscillatory,
    count_stable_fps, filter_stable_fps,
    count_oscillatory_fps


function jacobian_is_stable(jac::Matrix)
    all(real.(eigvals(jac)) .< 0)
end
function fixedpoint_is_stable(nullcline_params, fp)
    @assert all(abs.(derive_vector_fn(nullcline_params)(fp)) .< sqrt(eps())) derive_vector_fn(nullcline_params)(fp)
    jac_fn = derive_jacobian_fn(nullcline_params)
    jac = jac_fn(fp)
    jacobian_is_stable(jac)
end

function jacobian_is_oscillatory(jac::Matrix)
    any(imag.(eigvals(jac)) .!= 0)
end
function fixedpoint_is_oscillatory(nullcline_params::AbstractNullclineParams, fp)
    @assert all(abs.(derive_vector_fn(nullcline_params)(fp)) .< sqrt(eps())) derive_vector_fn(nullcline_params)(fp)
    jac_fn = derive_jacobian_fn(nullcline_params)
    jac = jac_fn(fp)
    jacobian_is_oscillatory(jac)
end

function fixedpoint_stability(nullcline_params, fp)
    jac_fn = derive_jacobian_fn(nullcline_params)
    jac = jac_fn(fp)
    if all(real.(eigvals(jac)) .< 0)
        return -1
    elseif all(real.(eigvals(jac)) .== 0)
        return 0
    else
        return 1
    end
end

function count_stable_fps(prototype_name::String, fixed_nt::NamedTuple, fp_arr::NamedAxisArray)
    prototype = get_prototype(prototype_name)
    map(TravelingWaveSimulationsPlotting.enumerate_nt(fp_arr)) do (nt, fps)
        params = get_nullcline_params(prototype(; fixed_nt..., nt...))
        count(TravelingWaveSimulationsPlotting.fixedpoint_is_stable.(Ref(params), fps))
    end
end

function count_oscillatory_fps(prototype_name::String, fixed_nt::NamedTuple, fp_arr::NamedAxisArray)
    prototype = get_prototype(prototype_name)
    map(TravelingWaveSimulationsPlotting.enumerate_nt(fp_arr)) do (nt, fps)
        params = get_nullcline_params(prototype(; fixed_nt..., nt...))
        count(TravelingWaveSimulationsPlotting.fixedpoint_is_oscillatory.(Ref(params), fps))
    end
end 

function filter_stable_fps(prototype_name::String, fixed_nt::NamedTuple, fp_arr::NamedAxisArray{NAMES}) where NAMES
    prototype = get_prototype(prototype_name)
    filtered_fp_arr_vec = map(TravelingWaveSimulationsPlotting.enumerate_nt(fp_arr)) do (nt, fps)
        params = get_nullcline_params(prototype(; fixed_nt..., nt...))
        filter_stable_fps(params, fps)
    end
    NamedAxisArray{NAMES}(reshape(filtered_fp_arr_vec, size(fp_arr)), axes_keys(fp_arr))
end
function filter_stable_fps(params::AbstractNullclineParams, fps::AbstractVector)
    filter(fp -> TravelingWaveSimulationsPlotting.fixedpoint_is_stable(params, fp), fps)
end
