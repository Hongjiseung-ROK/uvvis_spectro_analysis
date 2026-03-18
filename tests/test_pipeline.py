"""
test_pipeline.py — Smoke tests for the analysis pipeline.
Run: pytest tests/
"""

import sys
from pathlib import Path
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from uvvis.io import load_spectra
from uvvis.peak import detect_peaks, summarise_peaks
from uvvis.kinetics import fit_pseudo_first_order

DATA_CSV = ROOT / "data" / "processed" / "UV-Vis_results.csv"


@pytest.fixture(scope="module")
def data():
    return load_spectra(DATA_CSV)


class TestDataLoading:
    def test_wavelength_range(self, data):
        assert data.wl_rxn.min() < 260
        assert data.wl_rxn.max() > 490

    def test_time_points_present(self, data):
        assert set(data.times) == {0, 60, 120, 180, 240, 300, 360}

    def test_a0_peak_near_400nm(self, data):
        # 4-NP peak must be in 390–420 nm region
        mask = (data.wl_blank >= 390) & (data.wl_blank <= 420)
        assert data.abs_0[mask].max() > 1.5, "4-NP peak at t=0 should be > 1.5 AU"


class TestPeakDetection:
    def test_4np_decreases(self, data):
        results = summarise_peaks(data)
        abs_4np = [r.abs_4np for r in results if r.abs_4np is not None]
        # Absorbance should generally decrease
        assert abs_4np[0] > abs_4np[-1], "4-NP peak should decrease over time"

    def test_4ap_appears(self, data):
        results = summarise_peaks(data)
        # 4-AP peak should be detectable by t=120s
        r120 = next(r for r in results if r.time_s == 120)
        assert r120.abs_4ap is not None and r120.abs_4ap > 0.3

    def test_4np_wavelength_stable(self, data):
        results = summarise_peaks(data)
        wls = [r.wl_4np for r in results if r.wl_4np is not None and r.abs_4np > 0.1]
        # Peak position should not drift more than ±10 nm
        assert max(wls) - min(wls) < 10


class TestKinetics:
    def test_kapp_positive(self, data):
        result = fit_pseudo_first_order(data)
        assert result.k_app > 0

    def test_kapp_magnitude(self, data):
        result = fit_pseudo_first_order(data)
        # K_app should be in the range 0.005–0.030 s⁻¹ for this system
        assert 0.005 < result.k_app < 0.030, f"K_app={result.k_app:.5f} out of expected range"

    def test_r_squared(self, data):
        result = fit_pseudo_first_order(data, exclude_times=[300, 360])
        assert result.r_squared > 0.95, f"R²={result.r_squared:.4f} too low"

    def test_half_life_reasonable(self, data):
        result = fit_pseudo_first_order(data, exclude_times=[300, 360])
        # t½ should be 30–120 s for a typical Ag NP catalysed 4-NP reduction
        assert 30 < result.half_life < 120
