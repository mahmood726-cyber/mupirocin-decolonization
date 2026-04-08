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

## Pipeline Rules
See `GEMINI.md` for detailed research pipeline mandates.
