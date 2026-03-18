"""
peak.py — Peak detection and characterisation for UV-Vis spectra.

Experiment context
------------------
Two chemically distinct peaks are expected:
  1. 4-NP peak  (~400 nm): DECREASING over reaction time
                            4-nitrophenolate absorbs strongly at ~400 nm.
  2. 4-AP peak  (~300 nm): INCREASING then plateauing
                            4-aminophenol product absorbs at ~295-303 nm.

Peak behaviour summary (from experimental data):
  - t=0:   4-NP @ 399.95 nm, A=1.869  |  4-AP region unresolved
  - t=60s: 4-NP @ 396.95 nm, A=1.279  |  4-AP not yet distinct
  - t=120s:4-NP @ 398.96 nm, A=0.703  |  4-AP @ 302.93 nm, A=0.409
  - t=180s:4-NP @ 396.95 nm, A=0.309  |  4-AP @ 299.96 nm, A=0.492  ← max
  - t=240s:4-NP @ 398.03 nm, A=0.118  |  4-AP @ 301.99 nm, A=0.491
  - t=300s:4-NP absent                 |  4-AP @ 299.96 nm, A=0.450
  - t=360s:4-NP absent                 |  4-AP @ 296.04 nm, A=0.400  (slight decrease = product dilution / further rxn)
"""

from dataclasses import dataclass, field
from typing import Optional
import numpy as np
from scipy.signal import find_peaks


# Wavelength search windows (nm)
REGION_4NP = (390, 430)   # 4-nitrophenolate / 4-NP anion
REGION_4AP = (285, 320)   # 4-aminophenol product


@dataclass
class PeakResult:
    time_s: int
    # 4-NP peak
    wl_4np: Optional[float]
    abs_4np: Optional[float]
    # 4-AP peak
    wl_4ap: Optional[float]
    abs_4ap: Optional[float]
    # All detected local maxima
    all_peaks_wl: np.ndarray = field(default_factory=lambda: np.array([]))
    all_peaks_abs: np.ndarray = field(default_factory=lambda: np.array([]))

    def __repr__(self):
        return (
            f"PeakResult(t={self.time_s}s | "
            f"4-NP: {self.wl_4np:.1f} nm A={self.abs_4np:.4f} | "
            f"4-AP: {self.wl_4ap:.1f} nm A={self.abs_4ap:.4f})"
            if self.wl_4np and self.wl_4ap else
            f"PeakResult(t={self.time_s}s | partial data)"
        )


def _peak_in_region(wl: np.ndarray, ab: np.ndarray,
                    lo: float, hi: float) -> tuple[Optional[float], Optional[float]]:
    """Return (peak_wavelength, peak_absorbance) for the global max in a wavelength window."""
    mask = (wl >= lo) & (wl <= hi)
    if mask.sum() == 0:
        return None, None
    idx = np.argmax(ab[mask])
    return float(wl[mask][idx]), float(ab[mask][idx])


def detect_peaks(wl: np.ndarray, ab: np.ndarray, time_s: int) -> PeakResult:
    """
    Detect 4-NP and 4-AP peaks in a single spectrum.

    Parameters
    ----------
    wl     : wavelength array (nm), ascending
    ab     : absorbance array
    time_s : reaction time in seconds (used for labelling)
    """
    wl_4np, abs_4np = _peak_in_region(wl, ab, *REGION_4NP)
    wl_4ap, abs_4ap = _peak_in_region(wl, ab, *REGION_4AP)

    # All local maxima (for overview plotting)
    peak_idx, _ = find_peaks(ab, prominence=0.02, distance=8)

    return PeakResult(
        time_s=time_s,
        wl_4np=wl_4np,
        abs_4np=abs_4np,
        wl_4ap=wl_4ap,
        abs_4ap=abs_4ap,
        all_peaks_wl=wl[peak_idx],
        all_peaks_abs=ab[peak_idx],
    )


def summarise_peaks(spectra_data) -> list[PeakResult]:
    """Run detect_peaks over all time points in a SpectraData object."""
    results = []
    for t in spectra_data.times:
        wl = spectra_data.wavelength(t)
        ab = spectra_data.absorbance(t)
        results.append(detect_peaks(wl, ab, t))
    return results
