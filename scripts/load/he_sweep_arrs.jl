
using DrWatson, StaticArrays, AxisIndices, NamedDims

if !@isdefined(refresh_sweep_arrs)
    refresh_sweep_arrs = false
end

if !@isdefined(A_sweep_lower_bound)
    A_sweep_lower_bound = 0.1
end

if !@isdefined(A_sweep_upper_bound)
    A_sweep_upper_bound = 1.5
end


@time if !@isdefined(blocking_fp_arr) || !@isdefined(monotonic_fp_arr) || refresh_sweep_arrs
    
    nullcline_static_mods = (
        α=(0.4, 0.7), 
        firing_θI=0.2, blocking_θI=0.5, 
        n_lattice = 2,
        save_idxs=nothing, save_on=false, saveat=0.1
    ) 
    file, filename = produce_or_load(datadir(), Dict("lower" => A_sweep_lower_bound, "upper"=>A_sweep_upper_bound); 
        prefix = "blocking_fp_arr",
        suffix = "bson",
        force = refresh_sweep_arrs
    ) do c
        A_range = c["lower"]:0.1:c["upper"]
        nullcline_sweeping_mods = (Aee=A_range, Aei=A_range, Aie=A_range, Aii=A_range)
        blocking_nullcline_static_mods = nullcline_static_mods
        blocking_prototype_name = "full_dynamics_blocking"
        blocking_fp_arr = wcm_sweep_calculate_fixedpoints(
            blocking_prototype_name, 
            blocking_nullcline_static_mods,
            nullcline_sweeping_mods,
            ; 
            axis_length=100
        )
        return @dict(blocking_fp_arr, blocking_nullcline_static_mods, blocking_prototype_name)
    end
    @unpack blocking_fp_arr = file #, blocking_nullcline_static_mods, blocking_prototype_name = file
    file, filename = produce_or_load(datadir(), Dict("lower" => A_sweep_lower_bound, "upper"=>A_sweep_upper_bound); 
        prefix = "monotonic_fp_arr",
        suffix = "bson",
        force = refresh_sweep_arrs
    ) do c
        A_range = c["lower"]:0.1:c["upper"]
        nullcline_sweeping_mods = (Aee=A_range, Aei=A_range, Aie=A_range, Aii=A_range)
        monotonic_nullcline_static_mods = nullcline_static_mods
        monotonic_prototype_name = "full_dynamics_monotonic"
        monotonic_fp_arr = wcm_sweep_calculate_fixedpoints(
            monotonic_prototype_name, 
            monotonic_nullcline_static_mods,
            nullcline_sweeping_mods,
            ; 
            axis_length=100
        )
        return @dict(monotonic_fp_arr, monotonic_nullcline_static_mods, monotonic_prototype_name)
    end
    @unpack monotonic_fp_arr = file #, monotonic_nullcline_static_mods, monotonic_prototype_name = file
end

refresh_sweep_arrs = false