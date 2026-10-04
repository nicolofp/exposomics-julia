# exposomics-julia

Julia analyses of environmental exposures (the *exposome*), chemical mixtures
and omics layers (transcriptome, methylome, miRNA, metabolome), with an
emphasis on Bayesian models (BWQS, ARD, probabilistic PCA, variational
inference). Shared, dataset-agnostic code lives in the `ExposomicsJulia`
package (`src/`); each analysis is a Quarto document (`.qmd`, Julia engine).

This is exploratory code: methods and implementations may contain mistakes.

## Structure

```
exposomics-julia/
├── Project.toml                      # package ExposomicsJulia + analysis dependencies
├── src/
│   ├── ExposomicsJulia.jl            # module and exports
│   ├── transformations.jl            # rescaling, ECDF scores, missingness/imputation, z-scores
│   ├── stats.jl                      # OLS with inference, association scans, fit metrics, posterior summaries
│   ├── models.jl                     # Turing models: BWQS family, ARD, negative binomial, Gamma–Poisson, pPCA
│   └── plotting.jl                   # palette, group codes, embedding scatter
├── analyses/
│   ├── 00_data_ingestion/            # raw releases (R formats) → Arrow files in .data/ (R/knitr)
│   ├── 01_data_prep/                 # building analysis datasets
│   ├── 02_mixtures/                  # BWQS and hierarchical BWQS mixture models
│   ├── 03_bayesian_models/           # ARD, ADVI, RxInfer, negative binomial, Gamma–Poisson, GMM
│   ├── 04_dimensionality_reduction/  # autoencoders, t-SNE/UMAP, kernel PCA, pPCA, SVD/PCA
│   ├── 05_omics/                     # transcriptome, CpG, miRNA and metabolome analyses
│   └── 06_methods_sandbox/           # splines
├── .data/                            # input data (git-ignored)
└── results/                          # outputs (git-ignored)
```
