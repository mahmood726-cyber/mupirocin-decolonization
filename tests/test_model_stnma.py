"""Tests for the statistical core (src/model_stnma.py::run_mupirocin_stnma).

Covers the fail-closed input guards (fast, no MCMC) and one end-to-end fit on
the bundled fixture (runs a short MCMC chain, a few seconds).
"""
import copy
import importlib.util
import json
import os
import sys

import pytest

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_SRC = os.path.join(_ROOT, "src", "model_stnma.py")
_FIXTURE = os.path.join(_ROOT, "data", "mupirocin_synthesis_input.json")


def _load_model():
    sys.path.insert(0, os.path.dirname(_SRC))
    spec = importlib.util.spec_from_file_location("model_stnma_mod", _SRC)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_MODEL = _load_model()


def _fixture():
    with open(_FIXTURE, "r") as f:
        return json.load(f)


# --- Fast guard tests (no MCMC) ---------------------------------------------

def test_empty_decolonization_pool_raises():
    """No row with outcome == 'decolonization' must fail closed."""
    data = {"rcts": [
        {"arm": "A", "n": 100, "events": 50, "outcome": "MRSA BSI",
         "region": "North America"},
    ]}
    with pytest.raises(ValueError, match="decolonization"):
        _MODEL.run_mupirocin_stnma(data)


def test_zero_n_raises():
    """n <= 0 produces NaN/inf se; the function must reject it (F4 regression)."""
    data = {"rcts": [
        {"arm": "A", "n": 0, "events": 0, "outcome": "decolonization",
         "region": "North America"},
    ]}
    with pytest.raises(ValueError, match="n > 0"):
        _MODEL.run_mupirocin_stnma(data)


def test_events_exceed_n_raises():
    """events > n yields p > 1 (silently clamped, garbage se); must be rejected."""
    data = {"rcts": [
        {"arm": "A", "n": 5, "events": 10, "outcome": "decolonization",
         "region": "North America"},
    ]}
    with pytest.raises(ValueError, match="events"):
        _MODEL.run_mupirocin_stnma(data)


def test_valid_fixture_passes_guard():
    """The bundled fixture has only valid cells, so the guard must not trip."""
    data = _fixture()
    bad = [r for r in data["rcts"]
           if r["outcome"] == "decolonization"
           and (r["n"] <= 0 or r["events"] < 0 or r["events"] > r["n"])]
    assert bad == [], f"fixture unexpectedly has invalid decolonization cells: {bad}"


# --- End-to-end fit on the fixture (runs a short MCMC chain) -----------------

def test_full_fit_shape_and_bounds():
    results, summary, input_hash = _MODEL.run_mupirocin_stnma(_fixture())
    # Only North America has decolonization data, so exactly one region row.
    assert len(results) == 1
    row = results[0]
    assert row["region"] == "North America"
    # Probabilities must be ordered and strictly inside (0, 1).
    assert 0.0 < row["ci_low"] <= row["estimate"] <= row["ci_high"] < 1.0
    # Posterior for the Mupirocin+CHG arm should sit near the observed CLEAR
    # rate (850/1063 ~= 0.80); loose bound tolerant of MCMC noise.
    assert 0.6 < row["estimate"] < 0.9
    assert isinstance(input_hash, str) and len(input_hash) == 64


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
