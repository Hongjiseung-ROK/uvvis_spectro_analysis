#!/usr/bin/env python3
"""
run_analysis.py — Full UV-Vis analysis pipeline (CLI entry point)

Usage:
    python scripts/run_analysis.py

Outputs (written to output/):
    figures/fig1_spectra.png       — time-series UV-Vis spectra
    figures/fig2_kinetics.png      — pseudo first-order kinetics fit
    figures/fig3_peak_shift.png    — 4-NP decay / 4-AP formation
    reports/peak_summary.csv       — peak positions & absorbances per time point
    reports/kinetics_summary.csv   — K_app, R², t½, intercept
"""

import sys
from pathlib import Path

# Allow running from project root
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from uvvis.io import load_spectra
from uvvis.peak import summarise_peaks
from uvvis.kinetics import fit_pseudo_first_order
from uvvis.plot import plot_spectra, plot_kinetics, plot_peak_shift

import pandas as pd

# ── Paths ────────────────────────────────────────────────────────────────────
DATA_CSV    = ROOT / "data" / "processed" / "UV-Vis_results.csv"
FIG_DIR     = ROOT / "output" / "figures"
REPORT_DIR  = ROOT / "output" / "reports"

FIG_DIR.mkdir(parents=True, exist_ok=True)
REPORT_DIR.mkdir(parents=True, exist_ok=True)


def main():
    print("=" * 55)
    print("  UV-Vis Analysis Pipeline")
    print("  4-NP Catalytic Reduction by Ag NPs")
    print("=" * 55)

    # ── 1. Load data ──────────────────────────────────────────
    print("\n[1/4] Loading data ...")
    data = load_spectra(DATA_CSV)
    print(f"  Wavelength range (rxn): "
          f"{data.wl_rxn.min():.1f} – {data.wl_rxn.max():.1f} nm")
    print(f"  Time points: {data.times}")

    # ── 2. Peak detection ─────────────────────────────────────
    print("\n[2/4] Detecting peaks ...")
    peak_results = summarise_peaks(data)

    peak_rows = []
    for r in peak_results:
        print(f"  t={r.time_s:>3}s | "
              f"4-NP: {r.wl_4np or 'N/A':>8} nm  A={r.abs_4np or 'N/A':>8} | "
              f"4-AP: {r.wl_4ap or 'N/A':>8} nm  A={r.abs_4ap or 'N/A':>8}")
        peak_rows.append({
            "time_s":      r.time_s,
            "4NP_wl_nm":   r.wl_4np,
            "4NP_abs":     r.abs_4np,
            "4AP_wl_nm":   r.wl_4ap,
            "4AP_abs":     r.abs_4ap,
        })

    df_peaks = pd.DataFrame(peak_rows)
    peak_csv = REPORT_DIR / "peak_summary.csv"
    df_peaks.to_csv(peak_csv, index=False)
    print(f"  → saved: {peak_csv.relative_to(ROOT)}")

    # ── 3. Kinetics fit ───────────────────────────────────────
    print("\n[3/4] Fitting pseudo first-order kinetics ...")
    # Exclude t=300 and t=360 if 4-NP is near noise floor (optional)
    kinetics = fit_pseudo_first_order(
        data,
        peak_wl=399.95,
        exclude_times=[300, 360],   # near noise; remove to include all
    )
    print(kinetics.report())

    df_kinetics = pd.DataFrame([{
        "peak_wl_nm":     kinetics.peak_wavelength_nm,
        "A0":             kinetics.A0,
        "k_app_per_s":    kinetics.k_app,
        "k_app_x1e3_per_s": kinetics.k_app_milli,
        "intercept":      kinetics.intercept,
        "R_squared":      kinetics.r_squared,
        "half_life_s":    kinetics.half_life,
    }])
    kin_csv = REPORT_DIR / "kinetics_summary.csv"
    df_kinetics.to_csv(kin_csv, index=False)
    print(f"  → saved: {kin_csv.relative_to(ROOT)}")

    # ── 4. Plots ──────────────────────────────────────────────
    print("\n[4/4] Generating figures ...")

    fig1 = plot_spectra(data, output_path=FIG_DIR / "fig1_spectra.png")
    print(f"  → fig1_spectra.png")

    fig2 = plot_kinetics(kinetics, output_path=FIG_DIR / "fig2_kinetics.png")
    print(f"  → fig2_kinetics.png")

    fig3 = plot_peak_shift(peak_results, output_path=FIG_DIR / "fig3_peak_shift.png")
    print(f"  → fig3_peak_shift.png")

    print("\nDone. All outputs written to output/")
    print("=" * 55)


if __name__ == "__main__":
    main()
