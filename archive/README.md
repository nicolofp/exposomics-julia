# archive/

Scratch scripts kept verbatim because they contain *some* unique content but
are not exposomics analyses and were not converted to Quarto. They are **not**
covered by the top-level `Project.toml` (some need extra packages) and are
kept only for reference; nothing in `analyses/` depends on them.

| File | Origin | What it is | Extra packages |
|---|---|---|---|
| `met_classification/J_DuckDB.jl` | met_classification | DuckDB queries on a remote Arrow table and a local kallisto `abundance.tsv`, join, plus a few lines of MvNormal/MvTDist sampling | DuckDB |
| `met_classification/lathe_test.jl` | met_classification | Lathe.jl random-forest tutorial on `bwqs_tmp.csv` (Lathe API has changed since) | Lathe |
| `met_classification/decorrelation.jl` | met_classification | Eigen-decorrelation of an ill-conditioned matrix (scratch) | – |
| `met_classification/string_search_ts.jl` | met_classification | Brute-force best-matching subsequence search in simulated random-walk time series | – |

Delete this folder if none of it is worth keeping.
