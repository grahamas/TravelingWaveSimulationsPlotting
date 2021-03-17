
using StaticArrays, AxisIndices, NamedDims

if !@isdefined(refresh_sweep_arrs)
    refresh_sweep_arrs = false
end

@time if !@isdefined(blocking_fp_arr) || !@isdefined(monotonic_fp_arr) || refresh_sweep_arrs
    A_range = 0.1:0.1:1.5
    sweeping_mods = (Aee=A_range, Aei=A_range, Aie=A_range, Aii=A_range)
    static_mods = (
        α=(0.4, 0.7), 
        firing_θI=0.2, blocking_θI=0.5, 
        n_lattice = 2,
        save_idxs=nothing, save_on=true, saveat=0.1
    ) 
    file, filename = produce_or_load(datadir(), []; 
        prefix = "blocking_fp_arr",
        suffix = "bson",
        force = refresh_sweep_arrs
    ) do c
        blocking_fp_arr = wcm_sweep_calculate_fixedpoints(
            "full_dynamics_blocking", 
            static_mods,
            sweeping_mods,
            ; 
            axis_length=100
        )
        return @dict(blocking_fp_arr)
    end
    @unpack blocking_fp_arr = file
    file, filename = produce_or_load(datadir(), []; 
        prefix = "monotonic_fp_arr",
        suffix = "bson",
        force = refresh_sweep_arrs
    ) do c
        monotonic_fp_arr = wcm_sweep_calculate_fixedpoints(
            "full_dynamics_monotonic", 
            static_mods,
            sweeping_mods,
            ; 
            axis_length=100
        )
        return @dict(monotonic_fp_arr)
    end
    @unpack monotonic_fp_arr = file
end

refresh_sweep_arrs = false