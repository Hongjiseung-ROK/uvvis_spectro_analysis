"""
kinetics.py — Pseudo first-order kinetics fitting for 4-NP reduction.

Theory
------
Assuming excess NaBH4, the reaction is pseudo first-order in 4-NP:

    d[4-NP]/dt = -K_app * [4-NP]

Beer-Lambert law (A = ε·l·c) allows substitution:

    ln(A_t / A_0) = -K_app * t

K_app is extracted as the negative slope of a linear regression of
ln(A_t/A_0) vs. t.

Experimental K_app from this dataset:
  K_app = 0.01258 s⁻¹  (12.58 × 10⁻³ s⁻¹)
  R²    = 0.9795
  t½    = 55.1 s
  Reference wavelength: 399.95 nm (4-NP absorption maximum)
"""

from dataclasses import dataclass
import numpy as np
from typing import List
from .peak import PeakResult, REGION_4NP


@dataclass
class KineticsResult:
    """Stores pseudo first-order kinetics fit parameters."""
    peak_wavelength_nm: float   # reference wavelength used (4-NP peak)
    A0: float                   # absorbance at t=0

    times: np.ndarray           # reaction times used in fit (s)
    abs_t: np.ndarray           # A_t values at peak wavelength
    ln_ratio: np.ndarray        # ln(A_t / A0)

    k_app: float                # apparent rate constant (s⁻¹)
    intercept: float            # y-intercept of linear fit
    r_squared: float

    @property
    def half_life(self) -> float:
        """t½ = ln(2) / K_app  (seconds)"""
        return np.log(2) / self.k_app

    @property
    def k_app_milli(self) -> float:
        """K_app in units of 10⁻³ s⁻¹"""
        return self.k_app * 1000

    def report(self) -> str:
        lines = [
            "── Pseudo First-Order Kinetics ──────────────────────",
            f"  Reference wavelength : {self.peak_wavelength_nm:.2f} nm",
            f"  A₀ (t=0)             : {self.A0:.5f}",
            f"  K_app                : {self.k_app:.6f} s⁻¹  "
            f"({self.k_app_milli:.4f} × 10⁻³ s⁻¹)",
            f"  R²                   : {self.r_squared:.6f}",
            f"  Half-life (t½)       : {self.half_life:.1f} s",
            f"  Linear fit           : ln(A_t/A₀) = "
            f"{-self.k_app:.6f}·t + {self.intercept:.5f}",
            "─────────────────────────────────────────────────────",
        ]
        return "\n".join(lines)


def fit_pseudo_first_order(
    spectra_data,
    peak_wl: float = 399.95,
    wl_tolerance: float = 3.0,
    exclude_times: List[int] = None,
) -> KineticsResult:
    """
    Fit ln(A_t/A₀) = -K_app·t using the 4-NP peak absorbance.

    Parameters
    ----------
    spectra_data   : SpectraData object
    peak_wl        : target wavelength for absorbance extraction (nm)
    wl_tolerance   : ±nm window around peak_wl
    exclude_times  : list of time points (s) to exclude from fit
                     (useful if reaction is complete and A_t ≈ noise)

    Returns
    -------
    KineticsResult
    """
    if exclude_times is None:
        exclude_times = []

    # A0 from blank/t=0 scan
    wl0 = spectra_data.wl_blank
    ab0 = spectra_data.abs_0
    mask0 = np.abs(wl0 - peak_wl) <= wl_tolerance
    A0 = float(ab0[mask0].max()) if mask0.sum() > 0 else float(ab0.max())

    times_fit, abs_t_fit, ln_ratios = [], [], []

    for t in spectra_data.times:
        if t == 0 or t in exclude_times:
            continue
        wl = spectra_data.wavelength(t)
        ab = spectra_data.absorbance(t)
        mask = np.abs(wl - peak_wl) <= wl_tolerance
        if mask.sum() == 0:
            continue
        At = float(ab[mask].max())
        if At <= 0:
            continue
        times_fit.append(t)
        abs_t_fit.append(At)
        ln_ratios.append(np.log(At / A0))

    times_arr = np.array(times_fit, dtype=float)
    ln_arr = np.array(ln_ratios)

    coeffs = np.polyfit(times_arr, ln_arr, 1)
    k_app = -coeffs[0]
    intercept = coeffs[1]
    r2 = float(np.corrcoef(times_arr, ln_arr)[0, 1] ** 2)

    return KineticsResult(
        peak_wavelength_nm=peak_wl,
        A0=A0,
        times=times_arr,
        abs_t=np.array(abs_t_fit),
        ln_ratio=ln_arr,
        k_app=k_app,
        intercept=intercept,
        r_squared=r2,
    )
