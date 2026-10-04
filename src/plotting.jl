"Default colour palette (one colour per group)."
const PALETTE = [:coral, :yellow2, :lime, :turquoise2, :magenta, :dimgray]

"""
    group_codes(v) -> Vector{Int}

Map the distinct values of `v` (e.g. groups, sex, clusters) to
integer codes `1:k` in sorted order, for indexing [`PALETTE`](@ref).
"""
function group_codes(v)
    lv = sort(unique(v))
    lookup = Dict(x => i for (i, x) in enumerate(lv))
    return [lookup[x] for x in v]
end

"""
    scatter_embedding(E; groups = nothing, title = "", palette = PALETTE, kwargs...)

Scatter plot of a 2-D embedding. `E` can be `2 × N` (UMAP, Flux encoders)
or `N × 2` (t-SNE, SVD scores); the orientation is detected from the shape.
`groups` (optional) colours points by group.
"""
function scatter_embedding(E::AbstractMatrix; groups = nothing, title = "", palette = PALETTE, kwargs...)
    X = size(E, 1) == 2 && size(E, 2) != 2 ? permutedims(E) : E
    mc = groups === nothing ? palette[1] : palette[mod1.(group_codes(groups), length(palette))]
    return scatter(X[:, 1], X[:, 2]; mc = mc, title = title, label = "",
                   markersize = 3, markerstrokewidth = 0, kwargs...)
end
