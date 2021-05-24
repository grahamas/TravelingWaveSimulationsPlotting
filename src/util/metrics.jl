
function seizure_index(point::PT) where {PT<:AbstractVector}
    @assert all(1. .>= point .>= 0.)
    if sum(point) == 0.
        return 0.
    else
        (point[1] - point[2]) / (sum(point)) * maximum(point) # E-I/ E+I goes from -1 to 1
    end
end

function seizure_index(critical_points::AV) where {PT<:AbstractVector, AV<:AbstractVector{<:PT}}
    isempty(critical_points) && return 0
    maximum(seizure_index.(critical_points))
end

function contrast_metric(point::PT) where {PT<:AbstractVector}
    @assert all(1. .>= point .>= 0.)
    if sum(point) == 0.
        return 0.
    else
        (point[1] - point[2]) / (sum(point)) # E-I/ E+I goes from -1 to 1
    end
end

function contrast_metric(critical_points::AV) where {PT<:AbstractVector, AV<:AbstractVector{<:PT}}
    isempty(critical_points) && return 0
    maximum(contrast_metric.(critical_points))
end