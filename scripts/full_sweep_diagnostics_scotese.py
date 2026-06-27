"""
=============================================================================
full_sweep_diagnostics_scotese.py  —  Temporal-evolution plots (S&W workflow)
=============================================================================

Builds the headline time-dependent diagnostic figures from the per-slice
corrected NetCDFs in paths_scotese.CORRECTED_DIR.  Companion to
build_summary_stats_scotese.py:  stats produces the CSV/markdown, this
script produces the figures.

OUTPUTS  (in paths_scotese.OUTPUT_DIR)

  SW_full_sweep_diagnostics.png        4-panel:
      (1) 99th-percentile land elevation through time, input vs corrected
      (2) per-slice Δ RMS and Δ max
      (3) input continental hypsometry as a time-elevation heatmap
      (4) corrected continental hypsometry as a time-elevation heatmap

  SW_hypsometry_selected_ages.png      overlay of corrected hypsometric
      curves at selected S&W ages vs a schematic Earth-modern reference

  SW_metrics_by_era.png                per-era box/violin plots of the
      key diagnostics (bias before/after, RMS before/after, Δz RMS).
      Useful for paper figures.

INPUT
    paths_scotese.CORRECTED_DIR / "<age>Ma_corrected_SW.nc"
    paths_scotese.CORRECTED_DIR / "per_slice_stats_SW.csv"
        (produced by build_summary_stats_scotese.py)

USAGE
    cd <project>/scripts_Scotese
    python full_sweep_diagnostics_scotese.py

DEPENDENCIES
    numpy, pandas, netCDF4, matplotlib
=============================================================================
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np, pandas as pd, netCDF4 as nc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths_scotese import CORRECTED_DIR, OUTPUT_DIR, FIGURES_DIR

PER_SLICE_CSV = CORRECTED_DIR / "per_slice_stats_SW.csv"


def era_of(t_Ma: float) -> str:
    if t_Ma < 66:   return "Cenozoic"
    if t_Ma < 252:  return "Mesozoic"
    if t_Ma < 540:  return "Paleozoic"
    return "Proterozoic"


# Era boundaries used for shaded backgrounds (Ma)
ERA_BOUNDS = [(0, 66, "Cenozoic", "#fae6c8"),
              (66, 252, "Mesozoic", "#cfe6cf"),
              (252, 540, "Paleozoic", "#cfd9e6")]


def add_era_shading(ax, x_invert=True):
    """Shade era backgrounds on a time-axis plot.  Assumes x is age (Ma)."""
    for lo, hi, name, col in ERA_BOUNDS:
        ax.axvspan(lo, hi, color=col, alpha=0.5, zorder=0)
        ax.text((lo + hi) / 2, 1.02, name, transform=ax.get_xaxis_transform(),
                ha="center", va="bottom", fontsize=8, color="black")


def hyp_curve(arr, mask, n=200, land_only: bool = False):
    """Cumulative elevation distribution.
    If land_only=True, only cells with elevation > 0 are included — the
    cumulative axis runs from 0 (lowest land) to 1 (highest land)."""
    z = arr[mask]
    z = z[np.isfinite(z)]
    if land_only:
        z = z[z > 0]
    z = np.sort(z)
    if len(z) == 0: return np.full(n, np.nan)
    p = np.linspace(0, 1, len(z))
    pn = np.linspace(0, 1, n)
    return np.interp(pn, p, z)


def build_hypsometric_heatmap_data(land_only: bool = True):
    """Load all corrected NetCDFs, sort by age, and build hypsometry matrices.
    With land_only=True (default), each slice's hypsometry is computed over
    cells with elevation > 0 only — so the cumulative-area axis runs over
    the LAND portion of the continent at that age.
    Returns (ages, H_orig, H_corr)."""
    files = sorted(CORRECTED_DIR.glob("*Ma_corrected_SW.nc"),
                   key=lambda f: int(f.name.split("Ma")[0]))
    ages = []; H_orig = []; H_corr = []
    for f in files:
        try:
            with nc.Dataset(f) as d:
                t = int(round(float(d.target_age_Ma)))
                M = d.variables["M_orig"][:].astype(float)
                Mc = d.variables["M_corrected"][:].astype(float)
                cont = d.variables["continent_mask"][:].astype(bool)
        except Exception:
            continue
        ages.append(t)
        H_orig.append(hyp_curve(M, cont, 200, land_only=land_only))
        H_corr.append(hyp_curve(Mc, cont, 200, land_only=land_only))
    return np.array(ages), np.array(H_orig), np.array(H_corr)


# ---------------------------------------------------------------------------
# GMT-dem4-style colormap.  Matches the standard GMT dem4 CPT values
# (positive elevations only).  Used for the continental-hypsometry panels.
# ---------------------------------------------------------------------------
def make_gmt_dem4_cmap(top_m: float = 4000.0):
    """Matplotlib LinearSegmentedColormap that reproduces the GMT 'dem4' CPT.

    `top_m` is the elevation that maps to the cmap's upper end (matplotlib
    fractional position 1.0).  Pair the returned cmap with
    `Normalize(vmin=0, vmax=top_m)` so the absolute elevation->color
    mapping (green at sea level, brown at ~2 km, mauve near the top) is
    preserved across different value ranges.  Stops above `top_m` are
    clipped to fraction 1.0 (which displays as the colorbar's top
    triangle when `extend="max"`).
    """
    raw = [
        (   0.0, ( 95/255, 159/255,  95/255)),   # dark green at sea level
        ( 250.0, (140/255, 188/255,  90/255)),   # green
        ( 500.0, (220/255, 195/255,  65/255)),   # yellow-green
        (1000.0, (200/255, 170/255,  60/255)),   # olive
        (1500.0, (190/255, 130/255,  70/255)),   # tan
        (2000.0, (180/255,  95/255,  40/255)),   # red-brown
        (2500.0, (135/255,  85/255,  35/255)),   # brown
        (3000.0, (110/255,  90/255, 100/255)),   # greyish purple
        (3500.0, (155/255, 130/255, 170/255)),   # light mauve
        (4000.0, (200/255, 175/255, 210/255)),   # pale mauve
        (4500.0, (230/255, 215/255, 235/255)),   # very pale mauve
        (5000.0, (245/255, 240/255, 248/255)),   # near-white at top
    ]
    # Keep stops at or below top_m; clip any taller ones to fraction 1.0
    # so we don't lose the highest-band colour if top_m falls between
    # two named stops (and so the cmap has a defined value at frac=1).
    stops = [(min(e/top_m, 1.0), rgb) for e, rgb in raw if e <= top_m + 1e-6]
    if stops[-1][0] < 1.0:                            # ensure endpoint anchor
        stops.append((1.0, stops[-1][1]))
    cdict = {"red": [], "green": [], "blue": []}
    seen = set()
    for pos, (r, g, b) in stops:
        if pos in seen:                               # dedupe at frac=1 if needed
            continue
        seen.add(pos)
        cdict["red"].append((pos, r, r))
        cdict["green"].append((pos, g, g))
        cdict["blue"].append((pos, b, b))
    return mcolors.LinearSegmentedColormap("gmt_dem4", cdict, N=512)


def main():
    df = pd.read_csv(PER_SLICE_CSV).sort_values("t_Ma")
    ages_h, H_orig, H_corr = build_hypsometric_heatmap_data()

    # -------- Figure 1: 4-panel temporal evolution
    fig, axes = plt.subplots(4, 1, figsize=(13, 16))

    # (1) p99 land elevation
    ax = axes[0]
    ax.plot(df.t_Ma, df.land_orig_p99_m, label="S&W input", color="grey", lw=1.5)
    ax.plot(df.t_Ma, df.land_corr_p99_m, label="corrected", color="C3", lw=1.5)
    ax.set_ylabel("99th-percentile land elev (m)")
    ax.set_title("Top 1% of land elevation through time")
    ax.invert_xaxis(); ax.grid(True, alpha=0.3); ax.legend()
    add_era_shading(ax)

    # (2) correction magnitude
    ax = axes[1]
    ax.plot(df.t_Ma, df.delta_rms_m, label="Δ RMS over continent", color="C0", lw=1.5)
    ax.plot(df.t_Ma, df.delta_max_m, label="Δ max", color="C2", lw=1, alpha=0.7)
    ax.set_ylabel("correction magnitude (m)")
    ax.set_title("Correction magnitude through time")
    ax.invert_xaxis(); ax.grid(True, alpha=0.3); ax.legend()
    add_era_shading(ax)

    # (3, 4) — Continental hypsometry through time.  Land-only (z > 0);
    # 0-4000 m elevation scale; GMT dem4 colormap; values above the
    # 4 km top of the scale are extended with the colorbar's top triangle.
    HYPSO_TOP_M = 4000
    dem4 = make_gmt_dem4_cmap(top_m=HYPSO_TOP_M)
    norm = mcolors.Normalize(vmin=0, vmax=HYPSO_TOP_M)

    ax = axes[2]
    im = ax.imshow(H_orig.T, origin="lower",
                   extent=[ages_h.min(), ages_h.max(), 0, 1],
                   aspect="auto", cmap=dem4, norm=norm)
    ax.set_xlabel("age (Ma)"); ax.set_ylabel("cum. land area fraction")
    ax.set_title("Continental hypsometry through time — S&W input (land only)")
    ax.invert_xaxis()
    fig.colorbar(im, ax=ax, fraction=0.05, label="elevation (m)", extend="max")

    ax = axes[3]
    im = ax.imshow(H_corr.T, origin="lower",
                   extent=[ages_h.min(), ages_h.max(), 0, 1],
                   aspect="auto", cmap=dem4, norm=norm)
    ax.set_xlabel("age (Ma)"); ax.set_ylabel("cum. land area fraction")
    ax.set_title("Continental hypsometry through time — corrected (land only)")
    ax.invert_xaxis()
    fig.colorbar(im, ax=ax, fraction=0.05, label="elevation (m)", extend="max")

    plt.suptitle("Scotese & Wright assimilation — full-sweep diagnostics",
                 fontsize=14, y=0.995)
    plt.tight_layout()
    out = FIGURES_DIR / "SW_full_sweep_diagnostics.png"
    plt.savefig(out, dpi=140, bbox_inches="tight")
    print(f"wrote {out}")

    # -------- Figure 2: selected-age hypsometric curves (land only, 0-6000 m)
    fig, ax = plt.subplots(figsize=(10, 6))
    selected = [0, 50, 100, 200, 300, 400, 500, 540]
    cmap = plt.cm.plasma(np.linspace(0, 1, len(selected)))
    modern_p = np.linspace(0, 1, 200)
    for i, t in enumerate(selected):
        f = CORRECTED_DIR / f"{t}Ma_corrected_SW.nc"
        if not f.exists(): continue
        with nc.Dataset(f) as d:
            Mc = d.variables["M_corrected"][:].astype(float)
            cont = d.variables["continent_mask"][:].astype(bool)
        pc = hyp_curve(Mc, cont, 200, land_only=True)
        ax.plot(modern_p, pc, color=cmap[i], lw=1.6, label=f"{t} Ma")
    # Earth-modern reference (land only, 0-6000 m)
    ax.plot([0.0, 0.18, 0.42, 0.66, 0.82, 0.92, 0.97, 1.00],
            [0,    50,   300,  800,  1500, 2500, 3500, 5500],
            "--", color="black", label="Earth (modern, land only)", lw=1.2)
    ax.set_xlabel("cum. land area fraction"); ax.set_ylabel("elevation (m)")
    ax.set_title("Corrected continental hypsometry at selected ages — land only (Scotese & Wright)")
    ax.grid(True, alpha=0.3); ax.legend(loc="upper left", fontsize=8)
    ax.set_ylim(0, 6000)
    plt.tight_layout()
    out = FIGURES_DIR / "SW_hypsometry_selected_ages.png"
    plt.savefig(out, dpi=140, bbox_inches="tight")
    print(f"wrote {out}")

    # -------- Figure 3: per-era box plots
    df["era"] = df["t_Ma"].apply(era_of)
    metrics = [("bias_before_m", "Bias before (m)"),
               ("bias_after_m",  "Bias after (m)"),
               ("rms_before_m",  "RMS before (m)"),
               ("rms_after_m",   "RMS after (m)"),
               ("land_corr_p99_m", "Land p99 corrected (m)"),
               ("delta_rms_m",   "Δz RMS (m)")]
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    eras = ["Cenozoic", "Mesozoic", "Paleozoic"]
    for ax, (col, label) in zip(axes.flat, metrics):
        data = [df[df.era == e][col].dropna().values for e in eras]
        ax.boxplot(data, labels=eras, showfliers=False)
        ax.set_title(label, fontsize=11)
        ax.grid(True, alpha=0.3)
    plt.suptitle("S&W assimilation — diagnostics by geological era", fontsize=14)
    plt.tight_layout()
    out = FIGURES_DIR / "SW_metrics_by_era.png"
    plt.savefig(out, dpi=140, bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
