"""
diagnose_correction_distances.py
---------------------------------

For each Phanerozoic slice, measure the great-circle distance from
every continental cell with |Δz| > thresh to the nearest declustered
geochem sample.  Reports per-slice quantiles (50/75/90/95/99) of that
distance, plus the fraction of |Δz|>thresh area lying beyond various
radii (250 / 500 / 1000 / 2000 km).

This is a self-audit: the province-wise CDF rescaling can in principle
apply non-zero Δz at cells arbitrarily far from any sample (it only
needs the global province to have ≥ N_MIN_P samples).  The smoothed-
residual kernel has a 150-km Gaussian fall-off so is bounded.  The
diagnostic separates the two contributions where it can.

USAGE
    python diagnose_correction_distances.py             # all 109 slices
    python diagnose_correction_distances.py 50 100 200  # a subset

OUTPUTS
    data/corrected_Scotese/correction_distance_diag.csv
    outputs_Scotese/figures/correction_distance_diag.png
"""
from __future__ import annotations
import sys, argparse
from pathlib import Path
import numpy as np
import pandas as pd
import netCDF4 as nc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from paths_scotese import CORRECTED_DIR, OUTPUT_DIR
import assimilate_scotese as A
from sw_io import nearest_cell_index

R_EARTH_KM = 6371.0

DZ_THRESHOLDS = (100.0, 250.0, 500.0)      # m
RADIUS_BANDS = (250.0, 500.0, 1000.0, 2000.0)  # km

DIAG_CSV = CORRECTED_DIR / "correction_distance_diag.csv"
FIG_DIR  = OUTPUT_DIR / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)


def great_circle_km(lat1, lon1, lat2, lon2):
    """Vectorised — lat1/lon1 are scalars, lat2/lon2 arrays in degrees."""
    p1 = np.radians(lat1); p2 = np.radians(lat2)
    dlam = np.radians(lon2 - lon1)
    cosD = np.sin(p1) * np.sin(p2) + np.cos(p1) * np.cos(p2) * np.cos(dlam)
    cosD = np.clip(cosD, -1.0, 1.0)
    return R_EARTH_KM * np.arccos(cosD)


def nearest_sample_distance_km(cell_lats, cell_lons, sample_lats, sample_lons):
    """For each cell (lat, lon), return min great-circle distance to any sample.
    Loops over samples (typically 100-300), vectorises over cells (~5000-20000)."""
    if len(sample_lats) == 0:
        return np.full_like(cell_lats, np.nan)
    dmin = np.full_like(cell_lats, np.inf, dtype=float)
    for slat, slon in zip(sample_lats, sample_lons):
        d = great_circle_km(slat, slon, cell_lats, cell_lons)
        np.minimum(dmin, d, out=dmin)
    return dmin


def analyse_slice(t: int):
    """Returns (summary_row_dict, dist_array_km, dz_array_m) for one slice."""
    f = CORRECTED_DIR / f"{int(t)}Ma_corrected_SW.nc"
    with nc.Dataset(f) as d:
        lat   = np.asarray(d.variables["lat"][:])
        lon   = np.asarray(d.variables["lon"][:])
        M     = np.asarray(d.variables["M_orig"][:],         dtype=float)
        Mc    = np.asarray(d.variables["M_corrected"][:],    dtype=float)
        delta = np.asarray(d.variables["delta"][:],          dtype=float)
        cont  = np.asarray(d.variables["continent_mask"][:], dtype=bool)

    # Reconstruct the declustered-sample set used in the assimilation
    df = A.prepare_samples(A.get_geochem(), t)
    if df.empty:
        return {"t_Ma": int(t), "n_samples": 0,
                **{f"n_cells_dz_ge_{int(z)}m": 0 for z in DZ_THRESHOLDS}}, None, None
    iy = nearest_cell_index(lat, df["rlat"].values)
    ix = nearest_cell_index(lon, df["rlon"].values)
    df = df[cont[iy, ix]].reset_index(drop=True)
    if df.empty:
        return {"t_Ma": int(t), "n_samples": 0,
                **{f"n_cells_dz_ge_{int(z)}m": 0 for z in DZ_THRESHOLDS}}, None, None
    dec = A.decluster(df, lat, lon)
    sample_lats = dec["rlat"].values
    sample_lons = dec["rlon"].values

    # For every continental cell, compute distance to nearest declustered sample.
    # Restricting to continental cells: the polygon-vs-DEM ghost-cell issue
    # is muted by also requiring (M_orig >= 0) | (Mc >= 0).
    land = cont & ((M >= 0) | (Mc >= 0))
    cell_iy, cell_ix = np.where(land)
    cell_lats = lat[cell_iy]
    cell_lons = lon[cell_ix]
    dz_at_cells = delta[cell_iy, cell_ix]

    dist_km = nearest_sample_distance_km(cell_lats, cell_lons,
                                         sample_lats, sample_lons)

    row = {"t_Ma": int(t),
           "n_samples": int(len(dec)),
           "n_continental_cells": int(land.sum())}

    for z_thr in DZ_THRESHOLDS:
        mask = np.abs(dz_at_cells) >= z_thr
        n_corrected = int(mask.sum())
        row[f"n_cells_dz_ge_{int(z_thr)}m"] = n_corrected
        if n_corrected == 0:
            for q in (50, 75, 90, 95, 99):
                row[f"dist_p{q}_km_dz_ge_{int(z_thr)}m"] = np.nan
            for r in RADIUS_BANDS:
                row[f"frac_dz_ge_{int(z_thr)}m_beyond_{int(r)}km"] = np.nan
            continue
        d = dist_km[mask]
        for q in (50, 75, 90, 95, 99):
            row[f"dist_p{q}_km_dz_ge_{int(z_thr)}m"] = float(np.nanpercentile(d, q))
        for r in RADIUS_BANDS:
            row[f"frac_dz_ge_{int(z_thr)}m_beyond_{int(r)}km"] = float((d > r).mean())

    return row, dist_km, dz_at_cells


def main():
    p = argparse.ArgumentParser()
    p.add_argument("ages", nargs="*", type=int,
                   help="ages to diagnose (default: every available slice)")
    p.add_argument("--csv-only", action="store_true",
                   help="write per-slice CSV, skip the summary plot")
    args = p.parse_args()

    if args.ages:
        ages = sorted(args.ages)
    else:
        files = sorted(CORRECTED_DIR.glob("*Ma_corrected_SW.nc"),
                       key=lambda f: int(f.name.split("Ma")[0]))
        ages = [int(f.name.split("Ma")[0]) for f in files]
        # Skip the 0 Ma bypass slice (Δz ≡ 0)
        ages = [a for a in ages if a > 0]
    print(f"Diagnosing {len(ages)} slices …")

    rows = []
    for t in ages:
        row, _, _ = analyse_slice(t)
        rows.append(row)
        if all(np.isnan(row.get(f"dist_p50_km_dz_ge_{int(z)}m", np.nan))
               for z in DZ_THRESHOLDS):
            print(f"  {t:4d} Ma — n_samples={row['n_samples']:3d}  (no corrections)")
        else:
            print(f"  {t:4d} Ma — n_samples={row['n_samples']:3d}  "
                  f"|Δz|≥100m: n={row['n_cells_dz_ge_100m']:5d}, "
                  f"p50={row['dist_p50_km_dz_ge_100m']:5.0f} km, "
                  f"p90={row['dist_p90_km_dz_ge_100m']:5.0f} km   "
                  f"beyond 1000km: "
                  f"{100*row['frac_dz_ge_100m_beyond_1000km']:.0f} %")

    df = pd.DataFrame(rows).sort_values("t_Ma")
    df.to_csv(DIAG_CSV, index=False)
    print(f"\nwrote {DIAG_CSV}")

    if args.csv_only:
        return

    # Summary plot
    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

    ax = axes[0]
    for z_thr, colour in zip(DZ_THRESHOLDS, ["#2c7fb8", "#7fcdbb", "#edf8b1"]):
        ax.plot(df["t_Ma"], df[f"dist_p50_km_dz_ge_{int(z_thr)}m"],
                lw=1.4, color=colour, label=f"median, |Δz|≥{int(z_thr)} m")
        ax.fill_between(df["t_Ma"],
                        df[f"dist_p50_km_dz_ge_{int(z_thr)}m"],
                        df[f"dist_p90_km_dz_ge_{int(z_thr)}m"],
                        color=colour, alpha=0.18,
                        label=f"50–90th pct, |Δz|≥{int(z_thr)} m")
    ax.axhline(150, color="grey", lw=0.8, ls="--",
               label="150 km residual-kernel scale")
    ax.invert_xaxis()
    ax.set_ylabel("distance to nearest sample (km)")
    ax.set_title("Distance from corrected cells to the nearest declustered "
                 "geochem sample, per slice")
    ax.legend(fontsize=8, ncol=2, loc="upper left")
    ax.grid(alpha=0.3)

    ax = axes[1]
    for r, colour in zip(RADIUS_BANDS, ["#fdae6b", "#e6550d", "#a63603", "#8c2d04"]):
        ax.plot(df["t_Ma"], 100 * df[f"frac_dz_ge_100m_beyond_{int(r)}km"],
                lw=1.4, color=colour, label=f"|Δz|≥100 m beyond {int(r)} km")
    ax.invert_xaxis()
    ax.set_xlabel("age (Ma)")
    ax.set_ylabel("fraction of corrected area (%)")
    ax.set_title("Fraction of |Δz|≥100 m area lying farther than R from "
                 "the nearest sample")
    ax.legend(fontsize=8, ncol=2, loc="upper left")
    ax.grid(alpha=0.3)
    ax.set_ylim(0, 100)

    fig.tight_layout()
    out = FIG_DIR / "correction_distance_diag.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
