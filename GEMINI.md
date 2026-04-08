# GEMINI.md — mupirocin-decolonization Research Pipeline Rules

## Purpose
Ship a reproducible, OA-first Spatio-Temporal Network Meta-Analysis (ST-NMA) of **Mupirocin for S. aureus decolonization** as an **E156 micro-paper + GitHub repo + interactive HTML dashboard**.

## Session Workflow
Before you use any tool or make changes, briefly say what you're about to do, then do it. After each tool call, summarize what you found and what's next.

## Statistical Framework: ST-NMA
- **Methods:** Use Bayesian hierarchical models (PyMC) to synthesize RCT data (CLEAR, REDUCE MRSA, NICU trials) with regional burden and economic covariates.
- **Data Integration:** ClinicalTrials.gov, Figshare (CLEAR), Zenodo (NICU), PMC (SSTI decolonization).
- **Outcome:** Primary decolonization rate and MRSA recurrence.

## Non-negotiables
1. **OA-only**: no paywalls.
2. **No secrets**: redact before logs.
3. **Memory ≠ evidence**: certified claims must cite evidence locators + hashes.
4. **Fail-closed**: if validation incomplete, REJECT + reasons.
5. **Determinism**: fixed seeds, stable sorting, pinned versions.

## TruthCert (proof-carrying numbers)
- Every number must come from certified claims or be labeled UNCERTIFIED.
- Evidence locator + hash + transformation steps + validator outcomes.

## E156 Workbook Protection
- `C:\Users\user\OneDrive - NHS\Desktop\rewrite-workbook.txt`: Always keep `CURRENT BODY` up to date.
- **NEVER** modify or delete `YOUR REWRITE` sections.

## Quality Loop
- **Fix ALL issues in one pass.**
- **Test after EACH change.** Run full suite and report pass/fail count.

## SHIP Ritual
When ready to **SHIP**:
1. Run full test suite.
2. Demo on fixtures.
3. Perform TruthCert validation.
4. Generate E156 micro-paper (156 words max).
5. Deploy HTML dashboard to GitHub Pages.
6. Update master `INDEX.md` and workbook.

## Platform Defaults
- Python-first (PyMC); `python` not `python3`.
- Offline-first tests + fixtures.
