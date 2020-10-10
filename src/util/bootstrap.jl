

function bootstrap(fn::Function, naa::NamedAxisArray, kept_axes...; 
        min_prop, max_prop, n_samples)
    name_syms = _namedaxisarray_names(naa)
    subsampled_syms = Tuple(setdiff(name_syms, kept_axes))
    results_from_subsamples = map(fn, subsamples_over_axes(naa, subsampled_syms, min_prop, max_prop, n_samples))
end

struct NamedAxesSubsampleIdxs
    axes_dict
    min_N
    max_N
    n_samples
end
Base.length(nasi::NamedAxesSubsampleIdxs) = nasi.n_samples
Base.eltype(nasi::NamedAxesSubsampleIdxs) = typeof(first(nasi))

function get_subsample_sizes(subsampler, max_tries)
    all_dim_names = set(keys(subsampler.axes_dict))
    all_sample_sizes = subsampler.min_N:subsampler.max_N
    dim_names = copy(all_dim_names)
    sample_sizes = copy(all_sample_sizes)
    n_tries = 0
    subsample_sizes_dict = Dict()
    while true
        if n_tries > max_tries
            error("failed to find valid subsampling")
        end
        n_tries += 1
        while !isempty(dim_names)
            dim_name = sample(1, dim_names, replace=false)
            dim_len = length(subsampler.axes_dict[dim_name])
            dividing_sample_sizes = filter(1:dim_len) do dim_sample_size
                any(sample_sizes .% dim_sample_size .== 0)
            end 
            dim_sample_size = sample(1, dividing_sample_sizes)
            sample_sizes = sample_sizes[sample_sizes .% dim_sample_size .== 0] .÷ dim_sample_size
            subsample_sizes_dict[dim_name] = dim_sample_size
        end
        if prod(values(subsample_sizes_dict)) < subsampler.min_N
            subsample_sizes_dict = Dict()
            dim_names = copy(all_dim_names)
            sample_sizes = copy(all_sample_sizes)
        else
            @assert prod(values(subsample_sizes_dict)) <= subsampler.max_N
            return subsample_sizes_dict
        end
    end
end

# Return an index (dictionary, named idx) where the total size of the index
# is in the range min_N:max_N, subsampling the axes in axes_dict
# Note: subsample axes rather than coordinates so that you're comparing apples-to-apples:
#       the quantities we're bootstrapping are comparisons of averages. It only makes sense
#       to make those comparisons when the averages are over the same domain
function get_idxs(subsampler::NamedAxesSubsampleIdxs, max_tries=1000)
    subsample_sizes_dict = get_subsample_sizes(subsampler, max_tries)
    subsample_dict = Dict(name => sample(sample_size, subsampler.axes_dict[name]) 
        for (name, sample_size) in pairs(subsample_sizes_dict))
    return subsample_dict    
end
function iterate(subsampler::NamedAxesSubsampleIdxs, state=subsampler.n_samples)
    state == 0 ? nothing : (get_idxs(subsampler), state - 1)
end


function subsamples_over_axes(naa::NamedAxisArray, axis_names, 
        min_prop, max_prop, n_samples)
    axes_dict = Dict(name => axes_keys(naa)[dim(naa, axis_names)] for name in axis_names)
    N = prod(length.(values(axes_dict)))
    min_N, max_N = ceil(Int, min_prop * N), floor(Int, max_prop * N)
    @assert N >= min_N
    return NamedAxesSubsampleIdxs(axes_dict, min_N, max_N, n_samples)
end