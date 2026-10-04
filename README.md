# exposomics-julia

Julia analyses of environmental exposures (the *exposome*), chemical mixtures
and omics layers (transcriptome, methylome, miRNA, metabolome) in relation to
birth and child-health outcomes, with an emphasis on Bayesian models
(BWQS, ARD, probabilistic PCA, variational inference).

This repository merges and cleans up two older experiment repos:

* [`EnvBRAN`](https://github.com/nicolofp/EnvBRAN) – Jupyter notebooks (HELIX exposome, placenta metals, GEO methylation/miRNA, Turing/RxInfer models);
* [`met_classification`](https://github.com/nicolofp/met_classification) – loose `.jl` scripts (BWQS variants, ARD/negative-binomial regression, metabolomics, miRNA, JuMP).

Duplicated helpers and models now live once in a small package (`src/`), and
every analysis is a Quarto document (`.qmd`, Julia engine) without stored
outputs. As in the originals, this is exploratory code: methods and
implementations may contain mistakes – use at your own risk.

## Layout

```
exposomics-julia/
├── Project.toml                 # package ExposomicsJulia + all analysis deps (no Manifest)
├── src/
│   ├── ExposomicsJulia.jl       # module, exports
│   ├── data.jl                  # paths, loaders, rescale, ECDF scores, missingness/imputation, z-scores
│   ├── stats.jl                 # OLS with inference, exposure/omics-wide scans, fit metrics, posterior draws/summaries
│   ├── models.jl                # Turing models: BWQS family, ARD, NB, Gamma–Poisson, pPCA family
│   └── plotting.jl              # palette, group codes, embedding scatter
├── analyses/
│   ├── 01_data_prep/            # building analysis datasets (HELIX, MMIP)
│   ├── 02_mixtures/             # BWQS / hierarchical BWQS, metals and birth weight
│   ├── 03_bayesian_models/      # ARD, ADVI, RxInfer, negative binomial, Gamma–Poisson, GMM
│   ├── 04_dimensionality_reduction/  # autoencoders, t-SNE/UMAP, kernel PCA, pPCA, SVD/PCA of omics
│   ├── 05_omics/                # transcriptome, CpG, miRNA, metabolome association analyses
│   └── 06_methods_sandbox/      # splines, JuMP linear/integer programming
├── archive/                     # unconverted scratch scripts with some unique content
├── data/                        # git-ignored inputs (see data/README.md)
└── results/                     # git-ignored outputs (see results/README.md)
```

## Analyses

| Folder | Document | Topic |
|---|---|---|
| 01_data_prep | `helix_dataset_preparation.qmd` | HELIX pregnancy exposome + birth-weight z-score (cohort × gestational age × sex) |
| 01_data_prep | `mmip_data_casting.qmd` | MMIP methylation/exposure tables, long → wide |
| 02_mixtures | `helix_bwqs_metals.qmd` | BWQS of 9 prenatal metals on birth weight (flat vs hierarchical Dirichlet) |
| 02_mixtures | `bwqs_simulation.qmd` | BWQS variants (flat, hierarchical, logit-Dirichlet) on simulated mixtures; distributed chains |
| 02_mixtures | `hierarchical_bwqs.qmd` | Hierarchical BWQS with group intercepts/slopes |
| 02_mixtures | `placenta_metals_birthweight.qmd` | GSE71678 placenta metals vs birth weight: OLS and random forest |
| 03_bayesian_models | `ard_linear_regression.qmd` | ARD linear regression for exposome variable selection |
| 03_bayesian_models | `variational_inference.qmd` | ADVI for Bayesian linear regression |
| 03_bayesian_models | `rxinfer_exposure.qmd` | RxInfer message passing vs OLS |
| 03_bayesian_models | `negative_binomial_regression.qmd` | NB regression, ARD vs vague priors (simulation) |
| 03_bayesian_models | `gamma_poisson.qmd` | Gamma–Poisson count model (simulation) |
| 03_bayesian_models | `gaussian_mixture_turing.qmd` | Gaussian mixture with PG + HMC Gibbs |
| 04_dimensionality_reduction | `helix_autoencoders_tsne_umap.qmd` | Autoencoder, VAE, t-SNE, UMAP of exposures |
| 04_dimensionality_reduction | `helix_kernel_pca_umap_tsne.qmd` | Kernel PCA, UMAP, t-SNE of exposures/chemicals |
| 04_dimensionality_reduction | `probabilistic_pca.qmd` | pPCA / Dirichlet / ARD / supervised pPCA (simulated + ST001828) |
| 04_dimensionality_reduction | `svd_pca_omics.qmd` | SVD/PCA of serum metabolites, placental CpGs, miRNA (GSE48355) |
| 04_dimensionality_reduction | `low_rank_pca_metabolomics.qmd` | PCA/GLRM of ST001828 metabolites and GSE108497 miRNA |
| 05_omics | `helix_transcriptome_bmi.qmd` | Transcriptome-wide scan of child BMI |
| 05_omics | `helix_exposome_transcriptome_interactions.qmd` | Exposures, squares and interactions vs every transcript |
| 05_omics | `helix_urine_metabolome_elastic_net.qmd` | Exposure–metabolite correlations, lasso/elastic net |
| 05_omics | `gse71678_som_cpgs_metals.qmd` | SOM/t-SNE/UMAP of placental CpGs and metals |
| 05_omics | `gse154829_mirna_asthma.qmd` | miRNA and asthma: ARD logistic, Gamma–Poisson, SOM, DPMM |
| 05_omics | `gse59491_genexpr_preterm.qmd` | Maternal blood transcriptome, preterm birth, gestational age |
| 05_omics | `st001828_metabolomics_preterm.qmd` | Untargeted metabolomics: imputation, t/KS tests, random forest |
| 06_methods_sandbox | `splines.qmd` | B-spline fitting with BasicBSpline + Optim |
| 06_methods_sandbox | `linear_programming_jump.qmd` | Integer program, shortest path, min-cost flow with JuMP/HiGHS |

## Shared code (`src/`)

`using ExposomicsJulia` gives:

* **data** – `projectdir`, `datadir`, `resultsdir`, `read_arrow`, `load_helix_clean`, `HELIX_COVARIATES`, `helix_exposure_names`, `helix_table`, `rescale`, `ecdf_score(!)`, `missing_fraction`, `filter_by_missingness`, `impute_median`, `omics_matrix`, `add_zscore`;
* **stats** – `ols`, `ols_summary`, `exposure_scan`, `fit_metrics`, `draws`, `posterior_mean`, `chain_summary`, `credible_nonzero`;
* **models** – `bwqs`, `bwqs_adv`, `bwqs_soft`, `hbwqs`, `DirichletLogit`, `linear_regression`, `linear_regression_ard`, `logistic_ard`, `NegativeBinomial2`, `nb_regression_ard`, `negbin_regression`, `gamma_poisson`, `pPCA`, `pPCA_dir`, `pPCA_ARD`, `spPCA`, `spPCA_ARD`;
* **plotting** – `PALETTE`, `group_codes`, `scatter_embedding`.

## How to run

Requirements: Julia ≥ 1.10 (the environment was resolved and the package
loaded/smoke-tested with Julia 1.13.1 on Linux), Quarto ≥ 1.5 for rendering.

```bash
git clone https://github.com/nicolofp/exposomics-julia.git
cd exposomics-julia
julia --project=. -e 'using Pkg; Pkg.instantiate()'   # resolves a Manifest locally
```

Interactive use (from the repo root):

```julia
julia --project=.
julia> using ExposomicsJulia
julia> DT = load_helix_clean()     # needs data/helix/clean_dt.arrow
```

Put the input data under `data/` (layout and provenance in
[`data/README.md`](data/README.md)) or set `EXPOSOMICS_DATA=/path/to/data`.

**Quarto.** Every `.qmd` has

```yaml
engine: julia
julia:
  exeflags: ["--project=projects/exposomics-julia"]
```

because this repo is meant to be a git submodule at
`projects/exposomics-julia` of the main Quarto site, which renders from its
own root. There is deliberately no `_quarto.yml` here. Expensive or
data-dependent chunks (MCMC sampling, network training, batch scans, private
data) are marked `#| eval: false`; the rest needs the corresponding data in
`data/`. To render a page outside the site, run it from a directory where
`projects/exposomics-julia` points to this repo, or temporarily override the
flag, e.g. `quarto render analyses/06_methods_sandbox/splines.qmd -M julia.exeflags='["--project=."]'`
from the repo root.

## Notable changes vs. the original code

* **Ported to current package APIs** (the originals targeted 2022–23 versions):
  * Turing 0.49: `Turing.setadbackend`/`setrdcache` are gone – the AD backend is chosen in the sampler (`NUTS(; adtype = AutoReverseDiff(; compile = true))`); `sample` returns a FlexiChains chain, so `group(chain, :w)` / `chain.value[:, i, 1]` became `draws(chain, :w)`, `posterior_mean(chain, :w)` and `chain_summary(chain)`; `vi(m, ADVI(30, 1000))` became `vi(m, q_meanfield_gaussian, 1000)`; `Gibbs(PG(100, :k), HMC(0.05, 10, :μ, :w))` became `Gibbs(:k => PG(100), (:μ, :w) => HMC(0.05, 10))`.
  * Distributions: old `MvNormal(μ, σ)` calls (σ = standard deviation) were rewritten with explicit covariances (`MvNormal(μ, σ^2 * I)`, `Diagonal(s.^2)`), so the models are unchanged. Gamma–Poisson's `x[n, :] ~ Poisson(·)` is now `filldist(Poisson(·), D)`.
  * RxInfer 5: `datavar`/`inference`/`initmarginals` → data as model arguments, `infer`, `@initialization`.
  * Flux 0.16: implicit `Flux.params` + `train!(loss, ps, data, opt)` → `Flux.setup` + `train!(loss, model, data, state)`.
  * Convex 0.16: `silent_solver` → `silent`, `b.value` → `evaluate(b)`.
  * UMAP 0.2: `umap(X, k)` → `UMAP.fit(X, k).embedding`.
* **Evaluation.** Chunks that need data not in the repo, or that take minutes to sample, are marked `#| eval: false`; `helix_dataset_preparation`, `mmip_data_casting` and `placenta_metals_birthweight` set `execute: eval: false` for the whole document. With `data/helix/clean_dt.arrow` present, the remaining chunks were run end-to-end on Julia 1.13.1.
* **Bugs fixed (flagged here so results can be compared):**
  * `hbwqs` (from `bwqshd.jl`, `hbwqs2`): `(M * w) * βⱼ * τᵦ` summed the slopes of *all* groups for every observation; now each observation uses its own group's slope `βⱼ[idx] * τᵦ[idx]`.
  * `linear_regression_ard` (from `vs_bayesian.ipynb`): the coefficient prior mean was `ones(p)`; now `0` (`prior_mean = 1.0` reproduces the old model).
  * `logistic_ard`: dropped an unused `σ₂` parameter.
  * `NegativeBinomial2`: `p` is clamped to `[1e-4, 1 - 1e-10]` (the original only guarded `p > 0`; `p = 1` made NUTS fail at initialisation).
  * `ols`/`exposure_scan` return conventional **two-sided** p-values; the old loops used `cdf(TDist, -|t|)` (one-sided, i.e. half). Use `two_sided = false` to reproduce old numbers (done in `rxinfer_exposure.qmd`).
  * `GA_genexpr.jl` wrote into an undefined `RExp`; `EnvBRAN/dim_reduction.jl` coloured embeddings with vectors of the wrong length; `ppca_turing.ipynb`'s `pPCA_ARD` used `Normal(1.0, 0.0)`; `rca.jl` had hard-coded reshape sizes – all rewritten.
* **Column selection by name.** Positional selections that differed between the JLD2 and Arrow copies of the HELIX data (`17:73` vs `16:72`, `[1, 2, 6]`) use `helix_exposure_names(df)` and `HELIX_COVARIATES`; `load_helix_clean()` reads the Arrow file by default.
* **Dependencies trimmed.** `DPMMSubClusters` was removed from the environment (it pins JLD2 0.4/StatsBase 0.33 and held Turing, Flux and RxInfer back to versions that do not load on current Julia); its chunk in `gse154829_mirna_asthma.qmd` needs a separate environment. `LowRankModels` (incompatible with DataFrames ≥ 0.22) replaced by `MultivariateStats` PCA, the GLRM chunk kept as reference; `SyntheticDatasets` (needs Python) replaced by a 5-line circle generator; `Impute` replaced by `impute_median`; `GLPK` replaced by `HiGHS`; `Zygote`, `Memoization`, `AdvancedVI`, `FillArrays`, `MLDataUtils`, `MLBase`, `MLJ`, `NNlib`, `DuckDB`, `Lathe` dropped (unused or only in archived scripts).
* **Versions.** `[compat]` bounds come from a resolution done in Oct 2026 (Turing 0.49, RxInfer 5, Flux 0.16, DataFrames 1.8, JuMP 1.32, …). No Manifest is committed; `Pkg.instantiate()` resolves one locally. The old pinned `met_classification/Manifest.toml` was dropped.

## Old file → new location

### EnvBRAN

| Old file | New location / fate |
|---|---|
| `dataset_preparation.ipynb` | `analyses/01_data_prep/helix_dataset_preparation.qmd` |
| `b_mixture.ipynb` | `analyses/02_mixtures/helix_bwqs_metals.qmd` |
| `bw_placenta_metals.ipynb` | `analyses/02_mixtures/placenta_metals_birthweight.qmd` |
| `vs_bayesian.ipynb` | `analyses/03_bayesian_models/ard_linear_regression.qmd` (model → `src/models.jl`) |
| `variational_inference.ipynb` | `analyses/03_bayesian_models/variational_inference.qmd` (model → `src/models.jl`) |
| `rxinfer_exposure.md` | `analyses/03_bayesian_models/rxinfer_exposure.qmd` (the `.ipynb` only existed as an older checkpoint) |
| `test_gmm.ipynb` | `analyses/03_bayesian_models/gaussian_mixture_turing.qmd` |
| `dim_red_aut.ipynb` | `analyses/04_dimensionality_reduction/helix_autoencoders_tsne_umap.qmd` |
| `dim_reduction.jl` | `analyses/04_dimensionality_reduction/helix_kernel_pca_umap_tsne.qmd` |
| `ppca_turing.ipynb` | merged → `analyses/04_dimensionality_reduction/probabilistic_pca.qmd` |
| `eigen_svd.jl` | merged → `analyses/04_dimensionality_reduction/svd_pca_omics.qmd` |
| `Untitled.ipynb` | merged → `analyses/04_dimensionality_reduction/svd_pca_omics.qmd` (CpG SVD) |
| `As_miRNA.ipynb` (24 MB) | merged → `analyses/04_dimensionality_reduction/svd_pca_omics.qmd` (code only; outputs dropped) |
| `genexpr.jl` | `analyses/05_omics/helix_transcriptome_bmi.qmd` |
| `high_exposure.ipynb` | merged → `analyses/05_omics/helix_exposome_transcriptome_interactions.qmd` |
| `script_high_exposure.jl` | merged → same (batch section) |
| `umet_exposure.ipynb` | `analyses/05_omics/helix_urine_metabolome_elastic_net.qmd` |
| `gwas_test.ipynb` | merged → `analyses/05_omics/gse71678_som_cpgs_metals.qmd` |
| `j_arrow.ipynb` | merged → `analyses/05_omics/gse71678_som_cpgs_metals.qmd` |
| `splines.ipynb` | `analyses/06_methods_sandbox/splines.qmd` |
| `bayes_lib.jl` | merged → `src/models.jl` (`pPCA_ARD`, `spPCA_ARD`, `spPCA`, `bwqs`, `bwqs_adv`) |
| `README.md` | merged → this README |
| `*.md` exports (13: `b_mixture`, `bw_placenta_metals`, `dataset_preparation`, `dim_red_aut`, `gwas_test`, `high_exposure`, `j_arrow`, `ppca_turing`, `splines`, `test_gmm`, `umet_exposure`, `variational_inference`, `vs_bayesian` + rendered outputs) | dropped (rendered copies of the notebooks; code identical or older – checked cell by cell) |
| `*_files/` figure folders, `test_gmm.html` + `test_gmm_files/libs` | dropped (render artefacts) |
| `duckdb_test.ipynb` | dropped (2-cell stub with placeholder connection, no unique content) |
| `.ipynb_checkpoints/` (16 files) | dropped (older copies; `rxinfer_exposure-checkpoint` superseded by the `.md`) |
| `.gitignore` | replaced by the new `.gitignore` |
| `clead_dt.jld2`, `clean_dt.arrow` | not in git; copied locally to `data/helix/` |
| `mix_sim_chain.jld2`, `mix_adv_chain.jld2` | not copied (MCMC output; regenerated into `results/`) |
| `gexpr_results/`, `gexpr_results2/` (32,543 files) | not copied (batch output; regenerated into `results/gexpr_results/`) |

### met_classification

| Old file | New location / fate |
|---|---|
| `df_casting.jl` | `analyses/01_data_prep/mmip_data_casting.qmd` |
| `test_mixture.jl` | merged → `analyses/02_mixtures/bwqs_simulation.qmd` (models → `src/models.jl`) |
| `multichain_test.jl` | merged → `analyses/02_mixtures/bwqs_simulation.qmd` |
| `bayesian_stuff.jl` | split: BWQS part → `analyses/02_mixtures/bwqs_simulation.qmd`; NB part → `analyses/03_bayesian_models/negative_binomial_regression.qmd`; models → `src/models.jl` |
| `bwqshd.jl` | `analyses/02_mixtures/hierarchical_bwqs.qmd` (model `hbwqs` → `src/models.jl`) |
| `gp.jl` | `analyses/03_bayesian_models/gamma_poisson.qmd` |
| `ard_logistic_regression.jl` | merged → `src/models.jl` (`logistic_ard`) |
| `rca.jl` | merged → `analyses/04_dimensionality_reduction/probabilistic_pca.qmd` |
| `dim_reduction.jl` | merged → `analyses/04_dimensionality_reduction/low_rank_pca_metabolomics.qmd` |
| `clustering_cpgs.jl` | merged → `analyses/04_dimensionality_reduction/low_rank_pca_metabolomics.qmd` |
| `nb_regression.jl` | `analyses/05_omics/gse154829_mirna_asthma.qmd` |
| `GA_genexpr.jl` | `analyses/05_omics/gse59491_genexpr_preterm.qmd` |
| `00_data_cleaning.ipynb` | merged → `analyses/05_omics/st001828_metabolomics_preterm.qmd` (narrative) |
| `tmp_cleaning.jl` | merged → `analyses/05_omics/st001828_metabolomics_preterm.qmd` (newer, complete code) |
| `optimization_test.jl` | merged → `analyses/06_methods_sandbox/linear_programming_jump.qmd` |
| `network_minimal_flow.jl` | merged → `analyses/06_methods_sandbox/linear_programming_jump.qmd` |
| `simple_net.csv`, `simple_network_b.csv` | `analyses/06_methods_sandbox/network_data/` |
| `J_DuckDB.jl`, `lathe_test.jl`, `decorrelation.jl`, `string_search_ts.jl` | `archive/met_classification/` (see `archive/README.md`) |
| `Project.toml` | merged → `Project.toml` |
| `Manifest.toml` | dropped (environment re-resolved; no Manifest committed) |
| `README.md` | merged → this README |
| `tmp_julia.qmd` | dropped (Quarto "Plots Demo" template) |
| `gamma_poisson_1.txt` | dropped (empty) |
| `.ipynb_checkpoints/` | dropped |

### Deduplication into `src/`

| Helper / model | Copies replaced |
|---|---|
| `rescale` | 13 one-line definitions (two variants; the `eps()`-guarded one is canonical) |
| `ols` / `ols_summary` / `exposure_scan` | hand-written OLS + t-test loops in `bw_placenta_metals`, `high_exposure`, `script_high_exposure`, `genexpr.jl`, `GA_genexpr.jl`, `rxinfer_exposure` |
| `fit_metrics` | R²/adj-R²/MSE blocks (×2) in `bw_placenta_metals` |
| `ecdf_score!` | ECDF quantisation loops in `b_mixture`, `test_mixture`, `bayesian_stuff` |
| `add_zscore` | birth-weight z-score joins in `dataset_preparation`, `bw_placenta_metals` |
| `omics_matrix`, `impute_median`, `filter_by_missingness` | missing-filter + median-impute + transpose + rescale in `tmp_cleaning`, `dim_reduction` (met), `rca`, `clustering_cpgs` |
| `helix_table` | S3 TSV downloads in `umet_exposure`, `dim_reduction.jl` |
| `bwqs` | `bayes_lib.jl`, `test_mixture.jl`, `bayesian_stuff.jl` |
| `bwqs_adv` | `bwqs_adv` (`bayes_lib`) and `bwqs_new` (`test_mixture`, `bayesian_stuff`, `multichain_test`) – same model, different Gamma prior → keyword |
| `bwqs_soft`, `DirichletLogit`, `softmax` | `test_mixture.jl` (fixed α, ϕ) and `bayesian_stuff.jl` (learned) → keywords |
| `NegativeBinomial2`, `negbin_regression` | `NegativeBinomialRegression` and `negbinreg` (loop vs vectorised) |
| `gamma_poisson` | `gp.jl`, `nb_regression.jl` |
| `pPCA` | `ppca_turing.ipynb`, `rca.jl` |
| `chain_summary` | per-parameter quantile loops in `vs_bayesian`, `nb_regression`, `b_mixture` |
| `PALETTE`, `group_codes`, `scatter_embedding` | colour vectors and `recode` blocks in `dim_red_aut`, `j_arrow`, `As_miRNA` |
