# data/

**Git-ignored** (only this README is tracked). Put the input datasets here,
or point the environment variable `EXPOSOMICS_DATA` at the folder that holds
them (`datadir(...)` in `src/data.jl` resolves paths against it).

Expected layout, with where each dataset came from and the original local
path used in the old scripts:

| Folder / file | Dataset | Used by | Original location |
|---|---|---|---|
| `helix/clean_dt.arrow`, `helix/clead_dt.jld2` | Cleaned HELIX pregnancy exposome + birth-weight z-score (1004 × 72/73), built by `analyses/01_data_prep/helix_dataset_preparation.qmd` | 02, 03, 04, 05 | `EnvBRAN/clean_dt.arrow`, `EnvBRAN/clead_dt.jld2` (committed in EnvBRAN; local copies were placed here) |
| `helix/covariates.parquet`, `phenotype.parquet`, `codebook.parquet`, `genexpr_num.csv`, `genexpr_cov.parquet` | HELIX Exposome Data Challenge tables (ISGlobal), incl. transcriptome | 01, 05 | `C:/Users/nicol/Documents/Metabolomics/data/tables_format/` (not found on the PC in Oct 2026) |
| `helix/codebook.tsv`, `phenotype.tsv`, `covariates.tsv`, `met_urine.tsv` | Same HELIX tables as TSV, incl. urine metabolome | 04 (kernel PCA), 05 (elastic net) | `https://envbran.s3.us-east-2.amazonaws.com/<name>.tsv` (returns 403 as of Oct 2026); `helix_table()` reads the local copy first |
| `GSE71678/GSE71678_cov_original.arrow`, `GSE71678_cov_imputed.arrow`, `GSE71678_cpgs_original.arrow` | GEO GSE71678 – placenta DNA methylation, metals and birth outcomes (New Hampshire Birth Cohort) | 02, 04, 05 | `C:/Users/nicol/Documents/Datasets/GSE71678/Data/processed_data/` (also `Documents/GSE71678/...`) |
| `GSE48355/GSE48355_cov_original.arrow`, `GSE48355_miRNA.arrow`, `bw.csv` | GEO GSE48355 – placenta miRNA and arsenic; `bw.csv` = birth-weight reference (mean/SD by sex and week) | 04 | `C:/Users/nicol/Documents/Datasets/GSE48355/Data/` |
| `GSE154829/GSE154829_cov.arrow`, `GSE154829_miRNA.arrow` | GEO GSE154829 – miRNA and asthma | 05 | `C:/Users/nicol/Documents/Datasets/GSE154829/` |
| `GSE108497/GSE108497_miRNA_norm.arrow` | GEO GSE108497 – placental miRNA (normalised) | 04 | `C:/Users/nicol/Documents/Datasets/GSE108497/` |
| `GSE59491/GSE59491_cov.arrow`, `GSE59491_genexpr.arrow` | GEO GSE59491 – maternal blood gene expression, preterm birth | 05 | `C:/Users/nicol/Documents/Datasets/GSE59491/Data/` |
| `ST001828/covariates.arrow`, `bupos.arrow`, `buneg.arrow`, `btpos.arrow`, `btpos_code.arrow` | Metabolomics Workbench ST001828 – untargeted urine metabolomics, spontaneous preterm birth | 04, 05 | `C:/Users/nicol/Documents/Datasets/ST001828/` |
| `MMIP/epi.arrow`, `exposure.arrow`, `betas.arrow`, `mapping.arrow`, `codebook.arrow` | MMIP study (private) – epi data, chemical exposures, methylation | 01 | `G:/My Drive/Dati/MMIP/` |
| `bwqs_tmp.csv` | Small example dataset for BWQS (quantised urinary metals, covariates, cohort, `y`) | 02 | `https://envbran.s3.us-east-2.amazonaws.com/bwqs_tmp.csv` (403 now) or `C:/Users/nicol/Documents/bwqs_tmp.csv` |
| `serum/data_serum.json` | Serum metabolomics JSON | 04 | `https://provanik.s3.us-east-2.amazonaws.com/data_serum.json` (404 now) |

Only the two small HELIX files (`helix/clean_dt.arrow`, `helix/clead_dt.jld2`,
~1.4 MB) were copied here on the PC; the other sources were not found at
their original paths when the repo was created.
