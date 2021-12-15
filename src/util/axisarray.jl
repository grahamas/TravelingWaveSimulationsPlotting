_namedaxisarray_names(naa::NamedAxisArray{names}) where names = names
get_coordinates(aaa::AxisArray) = product(keys.(axes(aaa))...) |> collect

get_data(aa::AbstractArray) = aa
get_data(aa::AxisArray) = get_data(parent(aa))

get_axis_values(aa::NamedAxisArray, sym::Symbol) = keys.(axes(aa))[findfirst(AxisIndices.NamedDims.names(aa) .== sym)]
# should be:
# get_axis_values(aa::NamedAxisArray, sym::Symbol) = axes_keys(aa)[dim(aa, sym)]

nt_coord(naa::NamedAxisArray{X}, tup::Tuple) where X = NamedTuple{X}(tup)

macro ifsomething(ex)
    quote
        result = $(esc(ex))
        result === nothing && return nothing
        result
    end
end
struct EnumerateNAA{NAA,II}
    naa::NAA
    inner_iter::II
    EnumerateNAA(naa::NAA, ii::II) where {NAA,II} = new{NAA,II}(naa, ii)
end
EnumerateNAA(naa) = EnumerateNAA(naa, zip(Iterators.product(axes_keys(naa)...), naa))
Base.length(e::EnumerateNAA) = length(e.naa)

function enumerate_nt(naa::NamedAxisArray{X}) where X
    EnumerateNAA(naa)
end

function Base.iterate(enaa::EnumerateNAA, state)
    ((tup, val), state) = @ifsomething iterate(enaa.inner_iter, state)
    nt = nt_coord(enaa.naa, tup)
    return ((nt, val), state)
end
function Base.iterate(enaa::EnumerateNAA)
    ((tup, val), state) = @ifsomething iterate(enaa.inner_iter)
    nt = nt_coord(enaa.naa, tup)
    return ((nt, val), state)
end

function get_coordinate(naa::NamedAxisArray, cidx::CartesianIndex)
    get_coordinate(naa, Tuple(cidx))
end
function get_coordinate(naa::NamedAxisArray{NAMES}, tup::Tuple) where NAMES
    dim_coords = axes_keys(naa)
    NamedTuple{NAMES}(Tuple(dim[idx] for (dim, idx) ∈ zip(dim_coords, tup)))
end
function get_coordinate(naa::NamedAxisArray{NAMES}, idx::Int) where NAMES
    dim_coords = Iterators.product(axes_keys(naa)...) |> collect
    NamedTuple{NAMES}(dim_coords[idx])
end