# results/

**Git-ignored** (only this README is tracked). Analyses write cached MCMC
chains, tables and batch outputs here via `resultsdir(...)` (override the
location with `EXPOSOMICS_RESULTS`). Examples:

* `mix_sim_chain.jld2`, `mix_adv_chain.jld2` – BWQS chains of
  `analyses/02_mixtures/helix_bwqs_metals.qmd` (were committed in EnvBRAN);
* `gexpr_results/<gene>.txt` – per-transcript significant exposure terms of
  `analyses/05_omics/helix_exposome_transcriptome_interactions.qmd`
  (the ~32k `gexpr_results*/` files that lived in the EnvBRAN folder).
