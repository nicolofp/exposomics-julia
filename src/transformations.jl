"""
    rescale(A; dims = 1)

Centre and scale `A` along `dims` (z-score). Columns with zero variance are
divided by `eps()` instead of 0, so they become 0 rather than `NaN`.
"""
rescale(A; dims = 1) = (A .- mean(A; dims = dims)) ./ max.(std(A; dims = dims), eps())

"""
    ecdf_score(x; scale = 10)

Empirical-CDF score of a vector, multiplied by `scale` (`ecdf(x)(x) * scale`).
Used to turn exposures into continuous quantile scores for BWQS models.
"""
ecdf_score(x::AbstractVector; scale = 10) = ecdf(x)(x) .* scale

"""
    ecdf_score!(M; scale = 10)

Apply [`ecdf_score`](@ref) to every column of matrix `M` in place.
"""
function ecdf_score!(M::AbstractMatrix; scale = 10)
    for j in axes(M, 2)
        M[:, j] = ecdf_score(view(M, :, j); scale = scale)
    end
    return M
end

"""
    missing_fraction(M; dims = 2)

Number of missing entries along `dims` (`dims = 2` → per row, i.e. per feature
when features are rows as in the omics tables).
"""
missing_fraction(M; dims = 2) = vec(sum(ismissing.(M); dims = dims))

"""
    filter_by_missingness(df, id_col, value_cols; max_missing)

Keep the rows (features) of a features × samples table with at most
`max_missing` missing values across `value_cols`. Returns the filtered
`DataFrame` and the vector of kept feature ids.
"""
function filter_by_missingness(df::AbstractDataFrame, id_col, value_cols; max_missing::Integer)
    nmiss = missing_fraction(Matrix(df[:, value_cols]); dims = 2)
    keep = nmiss .<= max_missing
    return df[keep, :], string.(df[keep, id_col])
end

"""
    impute_median(M; dims = 2)

Replace missing values with the median of the non-missing values along
`dims` (`dims = 2` → row-wise, the convention used for features × samples
tables). Returns a `Float64` matrix.
"""
function impute_median(M::AbstractMatrix; dims::Integer = 2)
    out = Matrix{Float64}(undef, size(M))
    slices = dims == 2 ? axes(M, 1) : axes(M, 2)
    for i in slices
        v = dims == 2 ? view(M, i, :) : view(M, :, i)
        obs = collect(skipmissing(v))
        med = isempty(obs) ? NaN : median(obs)
        filled = [ismissing(x) ? med : Float64(x) for x in v]
        if dims == 2
            out[i, :] = filled
        else
            out[:, i] = filled
        end
    end
    return out
end

"""
    add_zscore(df, value_col, group_cols; name = Symbol(value_col, "_zscore"))

Add a column with the z-score of `value_col` computed within the strata
defined by `group_cols`. Intermediate `<value>_mean`/`<value>_std`
columns are kept.
"""
function add_zscore(df::DataFrame, value_col, group_cols; name = Symbol(value_col, "_zscore"))
    v = string(value_col)
    tab = combine(groupby(df, group_cols), v .=> [mean, std])
    joined = innerjoin(df, tab; on = group_cols)
    joined[!, name] = (joined[!, v] .- joined[!, v * "_mean"]) ./ joined[!, v * "_std"]
    return joined
end

"""
    omics_matrix(df, id_col, value_cols; max_missing = 0, standardize = true)

Turn a features × samples omics table (metabolites, miRNA, CpGs, …; one row
per feature, one column per sample) into an analysis matrix:

1. keep features with at most `max_missing` missing samples
   ([`filter_by_missingness`](@ref));
2. impute the remaining missing values with the feature median
   ([`impute_median`](@ref));
3. transpose to samples × features;
4. optionally standardise each feature ([`rescale`](@ref)).

Returns `(M, feature_ids)`.
"""
function omics_matrix(df::AbstractDataFrame, id_col, value_cols; max_missing::Integer = 0,
                      standardize::Bool = true)
    kept, ids = filter_by_missingness(df, id_col, value_cols; max_missing = max_missing)
    M = permutedims(impute_median(Matrix(kept[:, value_cols]); dims = 2))
    return (standardize ? rescale(M) : M), ids
end
