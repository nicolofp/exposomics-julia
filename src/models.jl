# Turing models.
#
# Choose the AD backend in the sampler, e.g.
# `NUTS(; adtype = AutoReverseDiff(; compile = true))`; `sample` returns a
# FlexiChains chain, see `draws`/`chain_summary` in stats.jl.

# ---------------------------------------------------------------------------
# Bayesian weighted quantile sum (BWQS) family
# ---------------------------------------------------------------------------

"""
    bwqs(cx, mx, y)

Bayesian weighted quantile sum regression:
`y ~ N(intercept + beta * (mx * w) + cx * delta, σ₂)` with
`w ~ Dirichlet(1, …, 1)`. `mx` holds quantile/ECDF scores of the mixture
components (see [`ecdf_score!`](@ref)), `cx` the covariates.
"""
@model function bwqs(cx, mx, y)
    σ₂ ~ Gamma(2.0, 2.0)
    intercept ~ Normal(0, 20)
    beta ~ Normal(0, 20)
    delta ~ filldist(Normal(0, 20), size(cx, 2))
    w ~ Dirichlet(size(mx, 2), 1)
    mu = intercept .+ beta * (mx * w) .+ cx * delta
    return y ~ MvNormal(mu, σ₂ * I)
end

"""
    bwqs_adv(cx, mx, y; alpha_prior = Gamma(1.0, 1.0))

BWQS with a hierarchical Dirichlet: `alpha_k ~ alpha_prior`, `w ~ Dirichlet(alpha)`.
"""
@model function bwqs_adv(cx, mx, y; alpha_prior = Gamma(1.0, 1.0))
    σ₂ ~ Gamma(2.0, 2.0)
    intercept ~ Normal(0, 20)
    beta ~ Normal(0, 20)
    delta ~ filldist(Normal(0, 20), size(cx, 2))
    alpha ~ filldist(alpha_prior, size(mx, 2))
    w ~ Dirichlet(alpha)
    mu = intercept .+ beta * (mx * w) .+ cx * delta
    return y ~ MvNormal(mu, σ₂ * I)
end

"Numerically naive softmax."
softmax(x) = exp.(x) ./ sum(exp.(x))

"""
    DirichletLogit(μ, ϕ)

Dirichlet parameterised by a location on the logit scale and a precision:
`Dirichlet(softmax(μ) * ϕ)`.
"""
DirichletLogit(μ, ϕ) = Dirichlet(softmax(μ) * ϕ)

"""
    bwqs_soft(cx, mx, y; alpha = nothing, phi = nothing)

BWQS with a logit-Dirichlet prior on the weights, `w ~ DirichletLogit(alpha, phi)`.
By default `alpha_k ~ Exponential(1)` and `phi ~ Gamma(2, 2)` are estimated.
Passing numbers fixes them, e.g. `alpha = fill(50.0, p), phi = 10.0`.
"""
@model function bwqs_soft(cx, mx, y; alpha = nothing, phi = nothing)
    σ₂ ~ Gamma(2.0, 2.0)
    intercept ~ Normal(0, 10)
    beta ~ Normal(0, 10)
    delta ~ filldist(Normal(0, 20), size(cx, 2))
    p = size(mx, 2)
    if alpha === nothing
        a ~ filldist(Exponential(1.0), p)
    else
        a = alpha
    end
    if phi === nothing
        ϕ ~ Gamma(2.0, 2.0)
    else
        ϕ = phi
    end
    w ~ DirichletLogit(a, ϕ)
    mu = intercept .+ beta * (mx * w) .+ cx * delta
    return y ~ MvNormal(mu, σ₂ * I)
end

"""
    hbwqs(M, X, idx, y)

Hierarchical BWQS (random intercept and random mixture slope per group
`idx ∈ 1:G`), non-centred slopes `βⱼ * τᵦ`. `M` = mixture scores, `X` =
covariates. Each observation uses its own group's slope.
"""
@model function hbwqs(M, X, idx, y; n_gr = length(unique(idx)), predictors = size(X, 2), nmix = size(M, 2))
    σ ~ Exponential(std(y))                               # residual SD
    τₐ ~ truncated(Cauchy(0, 2); lower = 0)               # SD of group intercepts
    δ ~ filldist(Normal(0, 20), predictors)               # fixed covariate effects
    αₖ ~ filldist(Gamma(2.0, 2.0), nmix)                  # Dirichlet concentration
    w ~ Dirichlet(αₖ)                                     # mixture weights
    τᵦ ~ filldist(truncated(Cauchy(0, 2); lower = 0), n_gr) # SD of group slopes
    αⱼ ~ filldist(Normal(0, τₐ), n_gr)                    # group intercepts
    βⱼ ~ filldist(Normal(0, 1), n_gr)                     # standardised group slopes
    ŷ = αⱼ[idx] .+ (M * w) .* (βⱼ[idx] .* τᵦ[idx]) .+ X * δ
    return y ~ MvNormal(ŷ, σ^2 * I)
end

# ---------------------------------------------------------------------------
# Regression models
# ---------------------------------------------------------------------------

"""
    linear_regression(x, y)

Bayesian linear regression with weakly-informative priors (used with ADVI in
the variational-inference analysis).
"""
@model function linear_regression(x, y)
    σ₂ ~ truncated(Normal(0, 100); lower = 0)
    intercept ~ Normal(0, sqrt(3))
    coefficients ~ MvNormal(zeros(size(x, 2)), 10.0 * I)
    mu = intercept .+ x * coefficients
    return y ~ MvNormal(mu, σ₂ * I)
end

"""
    linear_regression_ard(x, y; prior_mean = 0.0)

Linear regression with Automatic Relevance Determination:
`coefficients_j ~ N(prior_mean, 1/alpha_j)`, `alpha_j ~ Gamma(2, 2)`.
"""
@model function linear_regression_ard(x, y; prior_mean = 0.0)
    p = size(x, 2)
    σ ~ Gamma(1.0, 1.0)
    intercept ~ Normal(0, sqrt(3))
    alpha ~ filldist(Gamma(2.0, 2.0), p)
    coefficients ~ MvNormal(fill(prior_mean, p), Diagonal(1.0 ./ alpha))
    mu = intercept .+ x * coefficients
    return y ~ MvNormal(mu, σ^2 * I)
end

"""
    logistic_ard(x, y)

Bayesian logistic regression with ARD priors on the slopes
(`β₁_j ~ N(0, 1/α_j)`, `α_j ~ Gamma(2, 2)`).
"""
@model function logistic_ard(x, y)
    D = size(x, 2)
    α ~ filldist(Gamma(2.0, 2.0), D)
    β₀ ~ Normal(0, 20)
    β₁ ~ MvNormal(zeros(D), Diagonal(1.0 ./ α))
    v = StatsFuns.logistic.(β₀ .+ x * β₁)
    return y ~ arraydist(Bernoulli.(v))
end

"""
    NegativeBinomial2(μ, ϕ)

Negative binomial parameterised by mean `μ` and overdispersion `ϕ`
(variance `μ + μ²/ϕ`).
"""
function NegativeBinomial2(μ, ϕ)
    p = 1 / (1 + μ / ϕ)
    p = clamp(p, 1e-4, 1 - 1e-10) # numerical stability
    return NegativeBinomial(ϕ, p)
end

"""
    nb_regression_ard(x, y)

Negative-binomial regression (log link) with ARD priors on the slopes.
"""
@model function nb_regression_ard(x, y)
    D = size(x, 2)
    λ ~ Gamma(2.0, 2.0)
    α ~ filldist(Gamma(2.0, 2.0), D)
    β₀ ~ Normal(0, 20)
    β₁ ~ MvNormal(zeros(D), Diagonal(1.0 ./ α))
    mu = exp.(β₀ .+ x * β₁)
    return y ~ arraydist(LazyArray(@~ NegativeBinomial2.(mu, λ)))
end

"""
    negbin_regression(X, y)

Negative-binomial regression with vague Normal priors and `ϕ ~ InverseGamma(0.1, 0.1)`.
"""
@model function negbin_regression(X, y; predictors = size(X, 2))
    α ~ Normal(0, 10)
    β ~ filldist(Normal(0, 10), predictors)
    ϕ ~ InverseGamma(0.1, 0.1)
    return y ~ arraydist(LazyArray(@~ NegativeBinomial2.(exp.(α .+ X * β), ϕ)))
end

"""
    gamma_poisson(x)

Gamma–Poisson model for counts: row `n` of `x` is Poisson with a
row-specific background rate `br[n] ~ Gamma(a0, b0)`, with hyperpriors on
`a0`, `b0`.
"""
@model function gamma_poisson(x)
    N, D = size(x)
    a0 ~ Gamma(1, 1)
    b0 ~ Gamma(1, 1)
    br ~ filldist(Gamma(a0, b0), N)
    for n in 1:N
        x[n, :] ~ filldist(Poisson(br[n]), D)
    end
end

# ---------------------------------------------------------------------------
# Probabilistic PCA family
# ---------------------------------------------------------------------------

"""
    pPCA(X, k)

Probabilistic PCA (Turing tutorial form) for an `N × D` matrix `X` with `k`
latent dimensions: `X' ≈ W * Z .+ μ`.
"""
@model function pPCA(X::AbstractMatrix{<:Real}, k::Int)
    N, D = size(X)
    W ~ filldist(Normal(), D, k)
    Z ~ filldist(Normal(), k, N)
    μ ~ MvNormal(zeros(D), I)
    c_mean = W * Z .+ reshape(μ, D, 1)
    return X ~ arraydist([MvNormal(m, I) for m in eachcol(c_mean')])
end

"""
    pPCA_dir(X, k)

pPCA variant with Dirichlet (simplex) loadings: each column of `W` is a
probability vector over the `D` features.
"""
@model function pPCA_dir(X::AbstractMatrix{<:Real}, k::Int)
    N, D = size(X)
    W ~ filldist(Dirichlet(D, 1), k)
    Z ~ filldist(Normal(), k, N)
    μ ~ MvNormal(zeros(D), I)
    c_mean = W * Z .+ reshape(μ, D, 1)
    return X ~ arraydist([MvNormal(m, I) for m in eachcol(c_mean')])
end

"""
    pPCA_ARD(x, k)

Unsupervised pPCA with an ARD prior on the `k` latent components
(`alpha_k ~ Gamma(1, 1)`).
"""
@model function pPCA_ARD(x, k)
    N, D = size(x)
    K = Integer(k)
    z ~ filldist(Normal(), K, N)
    alpha ~ filldist(Gamma(1.0, 1.0), K)
    w ~ filldist(MvNormal(zeros(K), Diagonal((1.0 ./ alpha) .^ 2)), D)
    mu = (w' * z)'
    tau ~ Gamma(1.0, 1.0)
    for n in 1:N
        x[n, :] ~ MvNormal(mu[n, :], (1.0 / tau)^2 * I)
    end
end

"""
    spPCA(x, y, k)

Supervised pPCA: exposures `x` and outcome `y` share `k` latent factors.
"""
@model function spPCA(x, y, k)
    N, D = size(x)
    K = Integer(k)
    z ~ filldist(Normal(), K, N)
    w_x ~ filldist(Normal(), D, K)
    w_y ~ filldist(Normal(), K)
    mu_x = (w_x * z)'
    mu_y = (w_y' * z)'
    for n in 1:N
        x[n, :] ~ MvNormal(mu_x[n, :], I)
        y[n] ~ Normal(mu_y[n], 1.0)
    end
end

"""
    spPCA_ARD(x, y, m)

Supervised pPCA with ARD priors on the loadings of `x` and `y`.
"""
@model function spPCA_ARD(x, y, m)
    N, D = size(x)
    M = Integer(m)
    z ~ filldist(Normal(), M, N)
    alpha_x ~ filldist(Gamma(1.0, 1.0), M)
    alpha_y ~ Gamma(1.0, 1.0)
    w_x ~ filldist(MvNormal(zeros(M), Diagonal(alpha_x)), D)
    w_y ~ filldist(Normal(0.0, sqrt(alpha_y)), M)
    tau ~ Gamma(1.0, 1.0)
    for n in 1:N
        x[n, :] ~ MvNormal(w_x' * z[:, n], tau^2 * I)
        y[n] ~ Normal(w_y' * z[:, n], tau)
    end
end
