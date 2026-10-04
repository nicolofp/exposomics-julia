# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

"""
    projectdir(parts...)

Path relative to the root of this repository (works no matter which directory
Quarto/Julia is running from).
"""
projectdir(parts...) = joinpath(dirname(@__DIR__), parts...)

"""
    datadir(parts...)

Path inside the (git-ignored) data folder. Defaults to `<repo>/data`; set the
environment variable `EXPOSOMICS_DATA` to point somewhere else (e.g. a shared
drive holding the raw datasets).
"""
datadir(parts...) = joinpath(get(ENV, "EXPOSOMICS_DATA", projectdir("data")), parts...)

"""
    resultsdir(parts...)

Path inside the (git-ignored) results folder; created on demand. Override with
`EXPOSOMICS_RESULTS`.
"""
function resultsdir(parts...)
    root = get(ENV, "EXPOSOMICS_RESULTS", projectdir("results"))
    isdir(root) || mkpath(root)
    return joinpath(root, parts...)
end

# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

"""
    read_arrow(path) -> DataFrame

Read an Arrow file into a `DataFrame` (copying columns so the file handle can
be released).
"""
read_arrow(path::AbstractString) = DataFrame(Arrow.Table(path); copycols = true)

"""
    load_helix_clean(; format = :arrow)

Load the cleaned HELIX pregnancy-exposure dataset (1004 non-preterm births ×
covariates, birth-weight z-score and 57 pregnancy exposures) produced by
`analyses/01_data_prep/helix_dataset_preparation.qmd`.

* `format = :arrow` → `data/helix/clean_dt.arrow` (default; 72 columns)
* `format = :jld2`  → `data/helix/clead_dt.jld2` (73 columns: also has the
  categorical gestational-age class `e3_gac_Cat`, so every later column is
  shifted by one). Old JLD2 files may not load on recent Julia versions.

Select columns by name with [`HELIX_COVARIATES`](@ref) and
[`helix_exposure_names`](@ref) rather than by position.
"""
function load_helix_clean(; format::Symbol = :arrow)
    if format === :jld2
        return JLD2.load_object(datadir("helix", "clead_dt.jld2"))
    elseif format === :arrow
        return read_arrow(datadir("helix", "clean_dt.arrow"))
    else
        throw(ArgumentError("format must be :arrow or :jld2"))
    end
end

"Maternal covariates used as confounders in the HELIX birth-weight models (BMI, weight gain, age)."
const HELIX_COVARIATES = ["h_mbmi_None", "hs_wgtgain_None", "h_age_None"]

"""
    helix_exposure_names(df) -> Vector{String}

The 57 pregnancy exposures of the clean HELIX dataset (from
`h_abs_ratio_preg_Log` to `h_thm_preg_Log`), in table order. Replaces the
positional selections (`17:73` in the JLD2 file, `16:72` in the Arrow file,
`28:84` after joining the gene-expression covariates).
"""
function helix_exposure_names(df::AbstractDataFrame)
    n = names(df)
    i, j = findfirst(==("h_abs_ratio_preg_Log"), n), findfirst(==("h_thm_preg_Log"), n)
    (i === nothing || j === nothing) && throw(ArgumentError("HELIX exposure columns not found"))
    return n[i:j]
end

# ---------------------------------------------------------------------------
# Transformations
# ---------------------------------------------------------------------------

"""
    rescale(A; dims = 1)

Centre and scale `A` along `dims` (z-score). Columns with zero variance are
divided by `eps()` instead of 0, so they become 0 rather than `NaN`.
(Canonical version of the `rescale` one-liner defined in ~10 notebooks/scripts.)
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
tables). Returns a `Float64` matrix. Replaces `Impute.Substitute(; statistic = median)`.
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
    add_zscore(df, value_col, group_cols; name = :bw_zscore)

Add a column with the z-score of `value_col` computed within the strata
defined by `group_cols` (e.g. birth weight standardised by cohort,
gestational-age class and sex). Intermediate `<value>_mean`/`<value>_std`
columns are kept, matching the original notebooks.
"""
function add_zscore(df::DataFrame, value_col, group_cols; name = :bw_zscore)
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

Returns `(M, feature_ids)`. This is the preprocessing that was repeated in
`tmp_cleaning.jl`, `dim_reduction.jl`, `rca.jl` and `clustering_cpgs.jl`.
"""
function omics_matrix(df::AbstractDataFrame, id_col, value_cols; max_missing::Integer = 0,
                      standardize::Bool = true)
    kept, ids = filter_by_missingness(df, id_col, value_cols; max_missing = max_missing)
    M = permutedims(impute_median(Matrix(kept[:, value_cols]); dims = 2))
    return (standardize ? rescale(M) : M), ids
end

const HELIX_S3 = "https://envbran.s3.us-east-2.amazonaws.com"

"""
    helix_table(name) -> DataFrame

Read one of the HELIX tab-separated tables (`codebook`, `phenotype`,
`covariates`, `met_urine`, …) from `data/helix/<name>.tsv` if present,
otherwise try the original public bucket `$(HELIX_S3)/<name>.tsv`
(no longer publicly readable as of 2026 – keep a local copy).
"""
function helix_table(name::AbstractString)
    local_path = datadir("helix", name * ".tsv")
    src = isfile(local_path) ? local_path : HTTP.get("$(HELIX_S3)/$(name).tsv").body
    return CSV.read(src, DataFrame; delim = '\t')
end
