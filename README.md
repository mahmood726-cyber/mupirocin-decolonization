# mupirocin-decolonization: Spatio-Temporal NMA of S. aureus Decolonization

This repository implements a **Spatio-Temporal Network Meta-Analysis (ST-NMA)** to synthesize evidence for Mupirocin-based decolonization of *Staphylococcus aureus* (including MRSA).

## Project Scope
- **Data Integration:** Synthesizes RCT data from major decolonization trials (CLEAR, REDUCE MRSA, NICU phase 2) with regional antimicrobial resistance (AMR) burden and economic indicators.
- **Statistical Framework:** A Bayesian hierarchical model (PyMC) that models regional variations in decolonization efficacy and MRSA recurrence.
- **E156 Micro-Paper:** Includes a 7-sentence summary of findings with **TruthCert** proof-carrying numbers.

## Structure
- `src/`: Python scripts for data ingestion and Bayesian modeling.
- `data/`: Ingested (Open Access) data and fixed fixtures.
- `output/`: Model results, diagnostics, and TruthCert audit logs.
- `tests/`: Automated test suite for data validation and model stability.
- `docs/`: E156 micro-paper and project documentation.

## Deployment
Interactive dashboard hosted via GitHub Pages at `mahmood726-cyber.github.io/mupirocin-decolonization/`.

## Methods

The pooling model in `src/model_stnma.py` is a Bayesian arm-based hierarchical logistic NMA implemented in PyMC. Trial arms are mapped to intervention indices; observed decolonization probabilities are modelled on the logit scale with a per-intervention mean `mu_arm ~ Normal(0, 2)` and a per-region random effect `delta_region ~ Normal(0, tau_region)` with `tau_region ~ HalfNormal(0.5)`. Output is the posterior decolonization probability per region for the Mupirocin+CHG arm with 2.5%/97.5% credible bounds.

`data/mupirocin_synthesis_input.json` is the only input. The current trial set is **CLEAR**, **REDUCE MRSA**, and **NICU Phase 2**. The model filters on `outcome == "decolonization"`, so REDUCE MRSA — whose bundled outcome is MRSA BSI rather than decolonization — does not enter the decolonization pool; CLEAR and NICU Phase 2 do.

A SHA-256 `evidence_hash` is recorded per regional output. In the current implementation that hash is computed over the region name string only, not over the underlying trial cells; treat it as a deterministic tag, **not** as a cryptographic provenance certificate for the data.

## Limitations

- **Single-region trial coverage.** All three bundled trials are North America. Regional estimates for Europe, Asia, and Africa are produced from the model but are entirely prior-driven (`delta_region` has no likelihood contribution outside North America) — these should be read as "what the prior says about a region we have no evidence on", not as data-driven decolonization estimates.
- **MCMC settings are smoke-test, not production.** `pm.sample(200, tune=100, cores=1, chains=1)` is sufficient to validate the pipeline runs; it is **not** sufficient for trustworthy posterior summaries. With a single chain, Rhat cannot be computed and ESS will be far below the conventional ≥400-per-parameter threshold. Production runs should use at least 4 chains × 2000 post-warmup samples and verify `summary.r_hat < 1.01` and `summary.ess_bulk ≥ 400` before any reporting.
- **TruthCert hash is region-name only.** The `evidence_hash` field in each result row is `sha256(region_name)`. It does not capture trial cells, model code, or PyMC version. For real proof-carrying numbers, hash the (input JSON + model.py + library versions) tuple.
- **Outcome heterogeneity is hard-filtered, not modelled.** REDUCE MRSA's MRSA-BSI arm is dropped by the `outcome == "decolonization"` filter. Mixing surrogate / downstream outcomes (decolonization → infection) is not currently modelled; users wanting decolonization-as-mediator analyses need to extend the model.
- **`Mupirocin+CHG` is the silent reference.** If the input JSON lacks a 'Mupirocin+CHG' arm, `arm_map.get('Mupirocin+CHG', 0)` falls back to arm index 0 silently. Changing the input arms without updating this anchor will produce reference-arm confusion rather than an error.

## Conclusions

Use this repository as (a) a runnable scaffold for ST-NMA of decolonization studies and (b) an E156-protocol artifact demonstrating the pipeline. Do **not** quote regional estimates outside North America as evidence-driven, and do **not** treat the default MCMC settings as production. The next clear improvements are: (i) ingest at least one non-NA trial before claiming spatial variation, (ii) raise the sample/chain defaults and add an Rhat / ESS assertion gate, (iii) make the TruthCert hash cover the full (data, code, env) tuple.

## Pipeline Rules
See `GEMINI.md` for detailed research pipeline mandates.
