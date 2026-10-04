# ---------------------------------------------------------------------------
# Ordinary least squares with classical inference
# ---------------------------------------------------------------------------

"""
    ols(X, y; level = 0.95, two_sided = true)

Closed-form OLS fit (`β = X \\ y`) with standard errors, t statistics,
p-values and confidence intervals. `X` must already contain the intercept
column if one is wanted.

Returns a `NamedTuple` `(β, se, t, p, lower, upper, σ², dof, residuals)`.

!!! note
    The original notebooks computed `cdf(TDist(dof), -|t|)`, i.e. a
    **one-sided** p-value (half the usual two-sided value). The default here
    is the conventional two-sided p-value; pass `two_sided = false` to
    reproduce the numbers in the old notebooks.
"""
function ols(X::AbstractMatrix, y::AbstractVector; level = 0.95, two_sided::Bool = true)
    n, k = size(X)
    dof = n - k
    β = X \ y
    resid = y - X * β
    σ² = sum(abs2, resid) / dof
    Σ = σ² * inv(X' * X)
    se = sqrt.(diag(Σ))
    t = β ./ se
    tdist = TDist(dof)
    p = (two_sided ? 2 : 1) .* cdf.(tdist, -abs.(t))
    q = quantile(tdist, 1 - (1 - level) / 2)
    return (β = β, se = se, t = t, p = p, lower = β .- q .* se, upper = β .+ q .* se,
            σ² = σ², dof = dof, residuals = resid)
end

"""
    ols_summary(X, y, names; kwargs...) -> DataFrame

Coefficient table (`variable, beta, std_error, t_value, p_value, CI0025, CI0975`)
for [`ols`](@ref). `names` must have one entry per column of `X`.
"""
function ols_summary(X::AbstractMatrix, y::AbstractVector, names::AbstractVector; kwargs...)
    f = ols(X, y; kwargs...)
    return DataFrame(variable = string.(names), beta = f.β, std_error = f.se,
                     t_value = f.t, p_value = f.p, CI0025 = f.lower, CI0975 = f.upper)
end

"""
    exposure_scan(y, E, C; names = nothing, threaded = true, kwargs...) -> DataFrame

Exposome-/transcriptome-wide association scan: for every column `E[:, j]`
fit `y ~ 1 + E[:, j] + C` by OLS and keep the coefficient of `E[:, j]`.
`C` is a matrix of adjustment covariates (may have zero columns).

Adds a `bonferroni` column (`p_value < 0.05 / ncol(E)`).
This replaces the hand-written regression loops in `genexpr.jl`,
`high_exposure.ipynb`, `script_high_exposure.jl` and `GA_genexpr.jl`.
"""
function exposure_scan(y::AbstractVector, E::AbstractMatrix, C::AbstractMatrix;
                       names = nothing, threaded::Bool = true, alpha = 0.05, kwargs...)
    m = size(E, 2)
    T = promote_type(eltype(E), eltype(C), Float64)
    out = Matrix{Float64}(undef, m, 6)
    intercept = ones(T, size(E, 1))
    body = j -> begin
        X = hcat(intercept, E[:, j], C)
        f = ols(X, y; kwargs...)
        out[j, :] = [f.β[2], f.se[2], f.t[2], f.p[2], f.lower[2], f.upper[2]]
    end
    if threaded
        Threads.@threads for j in 1:m
            body(j)
        end
    else
        foreach(body, 1:m)
    end
    labels = names === nothing ? ["x$(j)" for j in 1:m] : string.(names)
    res = DataFrame(variable = labels, beta = out[:, 1], std_error = out[:, 2],
                    t_value = out[:, 3], p_value = out[:, 4],
                    CI0025 = out[:, 5], CI0975 = out[:, 6])
    res.bonferroni = res.p_value .< alpha / m
    return res
end

exposure_scan(y::AbstractVector, E::AbstractMatrix; kwargs...) =
    exposure_scan(y, E, zeros(size(E, 1), 0); kwargs...)

"""
    fit_metrics(y, ŷ, k) -> NamedTuple

R², adjusted R² and MSE for predictions `ŷ`, where `k` is the number of
estimated parameters (including the intercept).
"""
function fit_metrics(y::AbstractVector, ŷ::AbstractVector, k::Integer)
    n = length(y)
    ssres = sum(abs2, y .- ŷ)
    sstot = sum(abs2, y .- mean(y))
    r2 = 1 - ssres / sstot
    adjr2 = 1 - (ssres / (n - k)) / (sstot / (n - 1))
    return (R² = r2, adjR² = adjr2, MSE = ssres / n)
end

# ---------------------------------------------------------------------------
# Posterior summaries
# ---------------------------------------------------------------------------

"""
    credible_nonzero(lower, upper)

`true` where a credible/confidence interval excludes zero.
"""
credible_nonzero(lower, upper) = sign.(lower) .== sign.(upper)

"""
    draws(chain, name) -> Vector or Matrix

Posterior draws of parameter `name` (a `Symbol`, string or `VarName`) from a
Turing chain, pooling all chains. Scalars give a vector (one entry per draw);
vector/matrix-valued parameters give a `ndraws × length` matrix (column-major
order of the parameter). Works with the `FlexiChains` chains returned by
`sample` in Turing ≥ 0.40 (replaces `group(chain, :w)` / `chain[:, i, 1]`).
"""
function draws(chain, name)
    key = name isa AbstractString ? Symbol(name) : name
    v = vec(parent(chain[key]))
    return first(v) isa Number ? Float64.(v) : permutedims(reduce(hcat, [vec(x) for x in v]))
end

"""
    posterior_mean(chain, name)

Posterior mean of parameter `name`, keeping its shape (a scalar, vector or
matrix – e.g. the `k × N` latent scores `Z` of `pPCA`).
"""
posterior_mean(chain, name) = mean(vec(parent(chain[name isa AbstractString ? Symbol(name) : name])))

"""
    chain_summary(chain; quantiles = [0.025, 0.5, 0.975], pattern = nothing) -> DataFrame

Per-element posterior mean, sd and quantiles for every parameter of a Turing
chain (all chains pooled), plus a `nonzero` flag (outer quantiles exclude 0).
Vector/matrix parameters are expanded to `w[1]`, `Z[1,2]`, … . Restrict to
names containing `pattern` (e.g. `"coefficients"`, `"w["`). For convergence
diagnostics (ESS, R̂) use `Turing.FlexiChains.summarystats(chain)`.
"""
function chain_summary(chain; quantiles = [0.025, 0.5, 0.975], pattern = nothing)
    rows = NamedTuple[]
    for vn in Turing.FlexiChains.parameters(chain)
        vals = vec(parent(chain[vn]))
        base = string(vn)
        if first(vals) isa Number
            labels = [base]
            M = reshape(Float64.(vals), :, 1)
        else
            shape = size(first(vals))
            labels = [base * "[" * join(Tuple(I), ",") * "]" for I in CartesianIndices(shape)]
            M = permutedims(reduce(hcat, [vec(x) for x in vals]))
        end
        for (j, lab) in enumerate(labels)
            pattern === nothing || occursin(pattern, lab) || continue
            d = M[:, j]
            push!(rows, (; parameter = lab, mean = mean(d), sd = std(d),
                         (Symbol("q", q) => quantile(d, q) for q in quantiles)...))
        end
    end
    df = DataFrame(rows)
    qcols = [Symbol("q", q) for q in quantiles]
    df.nonzero = credible_nonzero(df[!, first(qcols)], df[!, last(qcols)])
    return df
end
