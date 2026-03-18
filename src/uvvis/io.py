"""
io.py — Data loading and preprocessing for UV-Vis CSV files.

CSV layout (two-row header, dual wavelength columns):
  Col A: wl_blank  | Col B: abs_0  (blank scan / t=0)
  Col C: wl_rxn    | Col D..I: abs_60 .. abs_360  (reaction series)
"""

import pandas as pd
import numpy as np
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict


REACTION_TIMES = [0, 60, 120, 180, 240, 300, 360]  # seconds


@dataclass
class SpectraData:
    """Container for parsed UV-Vis spectral data."""
    # t=0 / blank scan  (may span a slightly wider wavelength range)
    wl_blank: np.ndarray
    abs_0: np.ndarray

    # Reaction time series (common wavelength axis: wl_rxn)
    wl_rxn: np.ndarray
    spectra: Dict[int, np.ndarray]  # key = time in seconds

    @property
    def times(self):
        return sorted(self.spectra.keys())

    def absorbance(self, t: int) -> np.ndarray:
        """Return absorbance array for a given reaction time (s)."""
        if t == 0:
            return self.abs_0
        return self.spectra[t]

    def wavelength(self, t: int) -> np.ndarray:
        return self.wl_blank if t == 0 else self.wl_rxn


def load_spectra(csv_path: str | Path) -> SpectraData:
    """
    Load the UV-Vis CSV and return a SpectraData object.

    Parameters
    ----------
    csv_path : path to the processed CSV file

    Returns
    -------
    SpectraData
    """
    df = pd.read_csv(
        csv_path,
        header=None,
        skiprows=2,
        encoding="utf-8-sig",
    )
    df.columns = [
        "wl_blank", "abs_0",
        "wl_rxn",
        "abs_60", "abs_120", "abs_180",
        "abs_240", "abs_300", "abs_360",
    ]
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    # t=0 / blank: drop rows where wavelength is NaN
    df_blank = df[["wl_blank", "abs_0"]].dropna().sort_values("wl_blank").reset_index(drop=True)

    # Reaction series: rows where wl_rxn is present
    df_rxn = df[["wl_rxn", "abs_60", "abs_120", "abs_180",
                 "abs_240", "abs_300", "abs_360"]].dropna(
        subset=["wl_rxn"]
    ).sort_values("wl_rxn").reset_index(drop=True)

    spectra = {
        t: df_rxn[f"abs_{t}"].values
        for t in [60, 120, 180, 240, 300, 360]
    }

    return SpectraData(
        wl_blank=df_blank["wl_blank"].values,
        abs_0=df_blank["abs_0"].values,
        wl_rxn=df_rxn["wl_rxn"].values,
        spectra=spectra,
    )
