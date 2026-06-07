"""
=============================================================================
sw_io.py  —  Scotese & Wright NetCDF I/O helpers + age↔filename map
=============================================================================

The S&W PaleoDEM filenames encode the age (and sometimes a sub-index like
"Map80.5") instead of using a NetCDF time variable.  This module finds the
right file for a given age and provides a clean grid loader.

Public functions:
    available_ages()              -> sorted list[float] of all S&W ages
    file_for_age(age)             -> Path of the NetCDF for that age
    load_grid(age)                -> (M, lat, lon)
        M    : (n_lat, n_lon) float32, m  (variable is `z` in the file)
        lat  : 1-D, degrees, ascending (−90..+90 inclusive, 181 cells)
        lon  : 1-D, degrees, ascending (−180..+180 inclusive, 361 cells)

A small file-by-age index is built on first use and cached in memory.
"""
from __future__ import annotations
import re
from pathlib import Path
from functools import lru_cache
from typing import Iterable

import numpy as np
import netCDF4 as nc

from paths_scotese import DEM_DIR

# Filename pattern:  Map<index>[.<subindex>]_PALEOMAP_1deg_<period>_<age>Ma.nc
_AGE_RE = re.compile(r"_(\d+)Ma\.nc$")


@lru_cache(maxsize=1)
def _build_age_index() -> dict[int, Path]:
    """Return {age_Ma -> filepath} for all S&W PaleoDEM NetCDFs in DEM_DIR."""
    out: dict[int, Path] = {}
    for f in sorted(DEM_DIR.glob("*.nc")):
        m = _AGE_RE.search(f.name)
        if not m:
            continue
        age = int(m.group(1))
        # If multiple files match the same age (some have sub-indices like
        # Map80.5 vs Map80, with the same trailing "_460Ma"), prefer the one
        # without a sub-index — the .5-indexed maps are the auxiliary set.
        # Detect that by counting "." in the part before "_PALEOMAP"
        sub_indexed = "." in f.name.split("_PALEOMAP")[0]
        if age not in out or not sub_indexed:
            out[age] = f
    return out


def available_ages() -> list[int]:
    """Sorted list of ages (in Ma) for which an S&W PaleoDEM exists."""
    return sorted(_build_age_index().keys())


def file_for_age(age: float) -> Path:
    """Return the NetCDF Path for the given age (rounded to the nearest
    available age).  Raises FileNotFoundError if nothing within ±2.5 Myr."""
    idx = _build_age_index()
    target = int(round(age))
    if target in idx:
        return idx[target]
    # nearest within ±2.5 Myr
    nearest = min(idx.keys(), key=lambda a: abs(a - target))
    if abs(nearest - target) > 2.5:
        raise FileNotFoundError(
            f"No S&W PaleoDEM within ±2.5 Myr of {age} Ma "
            f"(nearest: {nearest} Ma).  Available: {available_ages()}"
        )
    return idx[nearest]


def load_grid(age: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read one S&W PaleoDEM at the given age.  Returns (M[lat,lon], lat, lon).

    The S&W grid is 181 lat × 361 lon (edge-aligned, both poles and both
    meridians included).  Cell *centres* therefore lie at integer degrees.
    """
    f = file_for_age(age)
    with nc.Dataset(f, "r") as d:
        z = d.variables["z"][:].astype(float)
        lat = d.variables["lat"][:].astype(float)
        lon = d.variables["lon"][:].astype(float)
    if hasattr(z, "filled"):
        z = z.filled(np.nan)
    # Ensure ascending lat (S&W files are ascending −90 … +90, but be safe)
    if lat[0] > lat[-1]:
        lat = lat[::-1]
        z = z[::-1, :]
    if lon[0] > lon[-1]:
        lon = lon[::-1]
        z = z[:, ::-1]
    return z, lat, lon


# ---------------------------------------------------------------------------
# Convenience: nearest-cell index lookup that works for either grid
# (cell-centred vs edge-aligned).  Used by the decluster step.
# ---------------------------------------------------------------------------
def nearest_cell_index(coord_array: np.ndarray, values: Iterable[float]) -> np.ndarray:
    """For each input value, return the index of the nearest cell centre.
    Works for arbitrary lat/lon ordering (must be monotonic ascending)."""
    v = np.asarray(values, dtype=float)
    # use searchsorted with edges at midpoints between consecutive cells
    midpoints = (coord_array[:-1] + coord_array[1:]) / 2.0
    idx = np.searchsorted(midpoints, v)
    return np.clip(idx, 0, len(coord_array) - 1)


if __name__ == "__main__":
    ages = available_ages()
    print(f"{len(ages)} S&W ages available: {ages[:5]} … {ages[-5:]}")
    print(f"  diffs (unique): {sorted(set([ages[i+1]-ages[i] for i in range(len(ages)-1)]))}")
    M, lat, lon = load_grid(50)
    print(f"50 Ma grid: shape={M.shape}, lat range {lat[0]}..{lat[-1]}, "
          f"lon range {lon[0]}..{lon[-1]}, z range {np.nanmin(M):.0f}..{np.nanmax(M):.0f}")
