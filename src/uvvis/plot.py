"""
plot.py — Publication-quality plotting for UV-Vis spectral analysis.

Figures produced:
  1. plot_spectra()      — Overlaid time-series UV-Vis spectra
  2. plot_kinetics()     — ln(A_t/A₀) vs. time with linear fit
  3. plot_peak_shift()   — 4-NP and 4-AP peak absorbance vs. time
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from pathlib import Path
from typing import Optional

from .io import SpectraData
from .kinetics import KineticsResult
from .peak import PeakResult, REGION_4NP, REGION_4AP


# ── Style constants ──────────────────────────────────────────────────────────
FONT_FAMILY = "DejaVu Sans"
LABEL_FONTSIZE = 12
TICK_FONTSIZE = 10
LEGEND_FONTSIZE = 9
DPI = 300

TIME_COLORS = {
    0:   "#2c3e50",   # near-black  (t=0, 4-NP dominant)
    60:  "#e74c3c",
    120: "#e67e22",
    180: "#f1c40f",
    240: "#2ecc71",
    300: "#1abc9c",
    360: "#3498db",   # blue        (4-AP dominant)
}


def _style_axes(ax, xlabel, ylabel, title=None):
    ax.set_xlabel(xlabel, fontsize=LABEL_FONTSIZE)
    ax.set_ylabel(ylabel, fontsize=LABEL_FONTSIZE)
    ax.tick_params(labelsize=TICK_FONTSIZE)
    if title:
        ax.set_title(title, fontsize=LABEL_FONTSIZE + 1, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


# ── Figure 1: Time-series spectra ────────────────────────────────────────────

def plot_spectra(
    data: SpectraData,
    output_path: Optional[Path] = None,
    show: bool = False,
) -> plt.Figure:
    """
    Overlaid UV-Vis spectra for all time points.

    Annotations:
      - Shaded region for 4-NP (~400 nm) in red
      - Shaded region for 4-AP (~300 nm) in blue
      - Arrow indicating spectral evolution direction
    """
    fig, ax = plt.subplots(figsize=(8, 5))

    for t in data.times:
        wl = data.wavelength(t)
        ab = data.absorbance(t)
        label = f"t = {t} s" if t > 0 else "t = 0 s (before rxn)"
        ax.plot(wl, ab, color=TIME_COLORS[t], linewidth=1.5,
                label=label, zorder=3)

    # Shade peak regions
    ax.axvspan(*REGION_4NP, alpha=0.08, color="#e74c3c",
               label="4-NP region (~400 nm)")
    ax.axvspan(*REGION_4AP, alpha=0.08, color="#3498db",
               label="4-AP region (~300 nm)")

    # Reference lines at known peak maxima
    ax.axvline(399.95, color="#e74c3c", linewidth=0.8, linestyle="--", alpha=0.6)
    ax.axvline(300.0, color="#3498db", linewidth=0.8, linestyle="--", alpha=0.6)

    # Annotations
    ax.annotate("4-NP\n(↓)", xy=(400, 1.5), fontsize=9,
                color="#e74c3c", ha="center", va="top")
    ax.annotate("4-AP\n(↑)", xy=(300, 0.55), fontsize=9,
                color="#3498db", ha="center", va="bottom")

    _style_axes(ax,
                xlabel="Wavelength (nm)",
                ylabel="Absorbance (AU)",
                title="UV-Vis Spectra: Catalytic Reduction of 4-NP by Ag NPs")

    ax.set_xlim(245, 510)
    ax.set_ylim(-0.05, 2.1)
    ax.legend(fontsize=LEGEND_FONTSIZE, loc="upper right",
              framealpha=0.7, ncol=2)
    fig.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=DPI, bbox_inches="tight")
    if show:
        plt.show()
    return fig


# ── Figure 2: Kinetics plot ──────────────────────────────────────────────────

def plot_kinetics(
    result: KineticsResult,
    output_path: Optional[Path] = None,
    show: bool = False,
) -> plt.Figure:
    """
    ln(A_t/A₀) vs. reaction time with linear fit overlay.
    """
    fig, ax = plt.subplots(figsize=(6, 5))

    # Data points
    ax.scatter(result.times, result.ln_ratio,
               color="#e74c3c", s=60, zorder=4, label="Experimental data")

    # Fitted line
    t_fit = np.linspace(result.times.min(), result.times.max(), 200)
    ln_fit = -result.k_app * t_fit + result.intercept
    ax.plot(t_fit, ln_fit, color="#2c3e50", linewidth=1.8,
            label=(f"Linear fit\n"
                   f"K$_{{app}}$ = {result.k_app*1000:.3f} × 10⁻³ s⁻¹\n"
                   f"R² = {result.r_squared:.4f}"))

    # Equation text
    eq_str = (f"ln(A$_t$/A$_0$) = {-result.k_app:.5f}·t + {result.intercept:.4f}")
    ax.text(0.05, 0.05, eq_str, transform=ax.transAxes,
            fontsize=9, color="#555555", va="bottom")

    _style_axes(ax,
                xlabel="Reaction Time (s)",
                ylabel="ln(A$_t$ / A$_0$)",
                title="Pseudo First-Order Kinetics (4-NP Reduction)")

    ax.legend(fontsize=LEGEND_FONTSIZE, loc="upper right")
    fig.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=DPI, bbox_inches="tight")
    if show:
        plt.show()
    return fig


# ── Figure 3: Peak absorbance vs. time ──────────────────────────────────────

def plot_peak_shift(
    peak_results: list[PeakResult],
    output_path: Optional[Path] = None,
    show: bool = False,
) -> plt.Figure:
    """
    Track the 4-NP and 4-AP peak absorbances over time on a dual-axis plot.
    Clearly shows the anti-correlated behaviour of reactant and product.
    """
    times_4np = [r.time_s for r in peak_results if r.abs_4np is not None]
    abs_4np   = [r.abs_4np for r in peak_results if r.abs_4np is not None]
    times_4ap = [r.time_s for r in peak_results if r.abs_4ap is not None]
    abs_4ap   = [r.abs_4ap for r in peak_results if r.abs_4ap is not None]

    fig, ax1 = plt.subplots(figsize=(7, 5))
    ax2 = ax1.twinx()

    ax1.plot(times_4np, abs_4np, "o-", color="#e74c3c",
             linewidth=2, markersize=7, label="4-NP (~400 nm, ↓)")
    ax2.plot(times_4ap, abs_4ap, "s--", color="#3498db",
             linewidth=2, markersize=7, label="4-AP (~300 nm, ↑)")

    ax1.set_xlabel("Reaction Time (s)", fontsize=LABEL_FONTSIZE)
    ax1.set_ylabel("Absorbance of 4-NP (AU)", color="#e74c3c",
                   fontsize=LABEL_FONTSIZE)
    ax2.set_ylabel("Absorbance of 4-AP (AU)", color="#3498db",
                   fontsize=LABEL_FONTSIZE)
    ax1.tick_params(axis="y", labelcolor="#e74c3c", labelsize=TICK_FONTSIZE)
    ax2.tick_params(axis="y", labelcolor="#3498db", labelsize=TICK_FONTSIZE)
    ax1.tick_params(axis="x", labelsize=TICK_FONTSIZE)
    ax1.set_title("Peak Absorbance vs. Time: 4-NP Decay & 4-AP Formation",
                  fontsize=LABEL_FONTSIZE + 1, fontweight="bold")

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2,
               fontsize=LEGEND_FONTSIZE, loc="center right")

    fig.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=DPI, bbox_inches="tight")
    if show:
        plt.show()
    return fig
