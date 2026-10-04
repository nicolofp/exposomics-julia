"""
    ExposomicsJulia

Shared helpers for the exposomics / omics analyses in `analyses/`.

The module collects code that used to be copy-pasted across the EnvBRAN
notebooks and the met_classification scripts:

* `data.jl`     – data paths, loading, standardisation, missing-data filtering,
                  quantile (ECDF) scoring of exposures, birth-weight z-scores.
* `stats.jl`    – fast OLS with inference, single-exposure (EWAS-style) scans,
                  fit metrics, posterior-chain summaries.
* `models.jl`   – canonical Turing models (BWQS family, ARD regressions,
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
import Arrow
import JLD2
import StatsFuns
import CSV
import HTTP

include("data.jl")
include("stats.jl")
include("models.jl")
include("plotting.jl")

# data.jl
export projectdir, datadir, resultsdir, read_arrow, load_helix_clean, HELIX_COVARIATES, helix_exposure_names
export rescale, ecdf_score, ecdf_score!, missing_fraction, filter_by_missingness, impute_median
export add_zscore, omics_matrix, helix_table
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
