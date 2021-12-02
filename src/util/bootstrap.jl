
struct Bootstrapped{T}
    results_from_subsampling::Vector{T}
    result_from_full::T
end

struct Estimated{T}
    mean::T
    band::Union{T,Tuple{T,T}}
end

function bootstrap(fn::Function, nda::NamedDimsArray{Names}, kept_axes...; 
        min_prop, max_prop, n_samples) where Names
    subsampled_syms = Tuple(setdiff(Names, kept_axes))
    results_from_subsamples = map(fn, subsamples_over_axes(nda, subsampled_syms, min_prop, max_prop, n_samples))
    result_from_full = fn(nda)
    return Bootstrapped(results_from_subsamples, result_from_full)
end

struct NamedAxesSubsampleIdxs
    array
    axes_dict
    min_N
    max_N
    n_samples
end
Base.length(nasi::NamedAxesSubsampleIdxs) = nasi.n_samples
Base.eltype(nasi::NamedAxesSubsampleIdxs) = typeof(first(nasi))

function get_subsample_sizes(subsampler, max_tries)
    all_dim_names = keys(subsampler.axes_dict) |> Set
    all_sample_sizes = subsampler.min_N:subsampler.max_N
    dim_names = copy(all_dim_names)
    sample_sizes = copy(all_sample_sizes)
    n_tries = 0
    subsample_sizes_dict = Dict()
    while true
        if n_tries > max_tries
            error("failed to find valid subsampling for $all_sample_sizes")
        end
        n_tries += 1
        while !isempty(dim_names)
            dim_name = sample(collect(dim_names), 1, replace=false) |> only
            pop!(dim_names, dim_name)
            dim_len = length(subsampler.axes_dict[dim_name])
            dividing_sample_sizes = filter(1:dim_len) do dim_sample_size
                any(sample_sizes .% dim_sample_size .== 0)
            end 
            dim_sample_size = sample(dividing_sample_sizes, 1) |> only
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
function get_idxs(subsampler::NamedAxesSubsampleIdxs, max_tries=500)
    subsample_sizes_dict = get_subsample_sizes(subsampler, max_tries)
    subsample_dict = Dict(name => sample(subsampler.axes_dict[name], sample_size)
        for (name, sample_size) in pairs(subsample_sizes_dict))
    return subsample_dict    
end
function Base.iterate(subsampler::NamedAxesSubsampleIdxs, state=subsampler.n_samples)
    state == 0 ? nothing : (getindex(subsampler.array; get_idxs(subsampler)...), state - 1)
end


function subsamples_over_axes
        nda::NamedDimsArray{Names}, 
        nda_dims::NamedTuple{Names}, 
        min_prop, max_prop, n_samples
    ) where Names
    axes_dict = Dict(name => axes_keys(naa)[dim(naa, name)] for name in axis_names)
    N = prod(length.(values(axes_dict)))
    @show axes_dict
    min_N, max_N = ceil(Int, min_prop * N), floor(Int, max_prop * N)
    @show min_N, max_N
    @info "calculating max_n_samples:"
    @info "$(binomial(N, max_N))"
    @assert N >= min_N
    return NamedAxesSubsampleIdxs(naa, axes_dict, min_N, max_N, n_samples)
end

import Makie: convert_arguments
convert_arguments(est::Estimated) = (est.estimate,)