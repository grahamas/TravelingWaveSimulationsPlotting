
function epilepsy_metric(point::PT) where {PT<:SVector{2}}
    @assert all(1. .>= point .>= 0.)
    if sum(point) == 0.
        return 0.
    else
        (point[1] - point[2]) / (sum(point)) * maximum(point) # E-I/ E+I goes from -1 to 1
    end
end

function epilepsy_metric(critical_points::AV) where {PT<:SVector{2}, AV<:AbstractVector{<:PT}}
    isempty(critical_points) && return 0
    maximum(epilepsy_metric.(critical_points))
end