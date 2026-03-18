"""
uvvis — UV-Vis spectroscopy analysis package
Experiment: Catalytic reduction of 4-NP → 4-AP by Ag NPs
Course: 2026-1 기기분석실험
"""

from .io import load_spectra
from .peak import detect_peaks, PeakResult
from .kinetics import fit_pseudo_first_order, KineticsResult
from .plot import plot_spectra, plot_kinetics, plot_peak_shift

__all__ = [
    "load_spectra",
    "detect_peaks", "PeakResult",
    "fit_pseudo_first_order", "KineticsResult",
    "plot_spectra", "plot_kinetics", "plot_peak_shift",
]
