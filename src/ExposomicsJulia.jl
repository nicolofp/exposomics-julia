"""
    ExposomicsJulia

Shared helpers for the exposomics / omics analyses in `analyses/`.

Contents:

* `transformations.jl` – standardisation, missing-data filtering and
                  imputation, quantile (ECDF) scoring, grouped z-scores.
* `stats.jl`    – fast OLS with inference, single-exposure (EWAS-style) scans,
                  fit metrics, posterior-chain summaries.
* `models.jl`   – Turing models (BWQS family, ARD regressions,
                  negative binomial, gamma–Poisson, probabilistic PCA).
* `plotting.jl` – small plotting helpers (colour palettes, embeddings).

Analyses load it with `using ExposomicsJulia` after activating this project
(`julia --project=.` from the repo root, or the Quarto `exeflags`).
"""
module ExposomicsJulia

using DataFrames
using Distributions
using LinearAlgebra
using Statistics
using StatsBase
using Turing
using LazyArrays
using Plots
import StatsFuns

include("transformations.jl")
include("stats.jl")
include("models.jl")
include("plotting.jl")

# transformations.jl
export rescale, ecdf_score, ecdf_score!, missing_fraction, filter_by_missingness, impute_median
export add_zscore, omics_matrix
# stats.jl
export ols, ols_summary, exposure_scan, fit_metrics, draws, posterior_mean, chain_summary, credible_nonzero
# models.jl
export bwqs, bwqs_adv, bwqs_soft, hbwqs, DirichletLogit, softmax
export linear_regression, linear_regression_ard, logistic_ard
export NegativeBinomial2, nb_regression_ard, negbin_regression, gamma_poisson
export pPCA, pPCA_dir, pPCA_ARD, spPCA, spPCA_ARD
# plotting.jl
export PALETTE, group_codes, scatter_embedding

end # module
