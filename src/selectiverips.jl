struct SelectiveRips{
    I<:Integer, T, A<:AbstractMatrix{T}
} <: Ripserer.AbstractRipsFiltration{T, Simplex}
    dist::A
    factor::T
    threshold::T
end

function SelectiveRips{I}(
    dist::AbstractMatrix{T}, factor;
    threshold=maximum(dist)
) where {I, T}
    issymmetric(dist) || throw(ArgumentError("`dist` must be a distance matrix"))
    return SelectiveRips{I, T, typeof(dist)}(dist, factor, T(threshold))
end
function SelectiveRips(dist, factor; threshold=maximum(dist))
    return SelectiveRips{Int}(dist, factor, threshold=threshold)
end

# interface implementation
Ripserer.n_vertices(rips::SelectiveRips) = size(rips.dist, 1)

@propagate_inbounds function Ripserer.dist(rips::SelectiveRips{T}, i, j) where T
    return ifelse(i == j, zero(T), rips.dist[i, j])
end

Ripserer.threshold(rips::SelectiveRips) = rips.threshold
Ripserer.birth(rips::SelectiveRips, i) = rips.dist[i, i]
Ripserer.simplex_type(::SelectiveRips{I, T}, dim) where {I, T} = Simplex{dim, T, I}

@propagate_inbounds function Ripserer.diam(
    flt::SelectiveRips, sx::Simplex{D}, us, v
) where D
    # Get the radius, same as in Rips.
    r = Ripserer.diam(sx)
    for u in us
        d = Ripserer.dist(flt, v, u)
        if ismissing(d) || d > threshold(flt)
            return missing
        else
            r = ifelse(r > d, r, d)
        end
    end
    if D ≥ 2
        n = 0 # number of edges with weight below factor * r
        vertices = (us..., v)
        for i in 1:length(vertices), j in i+1:length(vertices)
            d = Ripserer.dist(flt, vertices[j], vertices[i])
            n += d < min(r * flt.factor + (1 - flt.factor), r)
            if n ≥ D - 1
                return r
            end
        end
        return missing
    else
        return r
    end
end
