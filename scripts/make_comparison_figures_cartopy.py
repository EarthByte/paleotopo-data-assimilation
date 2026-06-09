"""
=============================================================================
make_comparison_figures_cartopy.py  —  Figures 5a/b/c of the Earth-Science Reviews paper
=============================================================================

Produces three 2-column × 3-row comparison figures, one per pair of ages:
    Fig05a: 500 Ma  |  400 Ma
    Fig05b: 300 Ma  |  200 Ma
    Fig05c: 100 Ma  |   50 Ma

Layout (per figure):
    row 1 — Scotese & Wright 2018 original elevation
    row 2 — geochemically-corrected elevation + reconstructed declustered
            sample overlay (filled circles, coloured by observed
            elevation z_obs on the same terrain cmap as the background)
    row 3 — Δz (corrected − original) + reconstructed declustered
            sample overlay (small open black circles marking sample
            positions only — observed elevation is not directly
            comparable to Δz, so colour-coding by it would mislead)
    columns — the two ages

The samples shown are the per-slice declustered observations that
actually entered the assimilation — i.e. one weighted-median value
per (1° cell × Tecto_Prov) bin, reconstructed to the target age via
the S&W rotation model and dropped if they land offshore or on a
flooded-land cell.  Reader can therefore see directly where the
corrections came from.

Rendered in matplotlib + cartopy Robinson projection (the closest
pseudocylindrical to Winkel-Tripel available without GMT).  Fixed cpt:
    elevation  : terrain, −4000 … +4000 m, end-extension arrows
    Δz         : RdBu_r,  ±2000 m, end-extension arrows

OUTPUT  (paper/Scotese/)
    Fig05a_SW_comparison_500-400Ma.png
    Fig05b_SW_comparison_300-200Ma.png
    Fig05c_SW_comparison_100-50Ma.png

(A pyGMT version for publication-grade Winkel-Tripel rendering is in
`make_comparison_figures.py`; run that locally where GMT is installed.)
=============================================================================
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import netCDF4 as nc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from paths_scotese import CORRECTED_DIR
import assimilate_scotese as A
from sw_io import nearest_cell_index
PROJ_ROOT = HERE.parent
# Paper-numbered output goes to paper/Scotese/.  The pyGMT version of the
# same figure (make_comparison_figures.py) writes its own Fig05a/b/c PNGs
# straight to paper/Scotese/ via the build script — keeping these here too
# means either renderer leaves the figures in their final home.
FIG_DIR = PROJ_ROOT / "paper" / "Scotese"
FIG_DIR.mkdir(parents=True, exist_ok=True)

ELEV_NORM = mcolors.TwoSlopeNorm(vmin=-4000, vcenter=0, vmax=4000)
DELTA_NORM = mcolors.TwoSlopeNorm(vmin=-2000, vcenter=0, vmax=2000)

PAIRS = [(500, 400, "Fig05a"),
         (300, 200, "Fig05b"),
         (100,  50, "Fig05c")]


def load(t):
    f = CORRECTED_DIR / f"{int(t)}Ma_corrected_SW.nc"
    with nc.Dataset(f) as d:
        lat   = np.asarray(d.variables["lat"][:])
        lon   = np.asarray(d.variables["lon"][:])
        M     = np.asarray(d.variables["M_orig"][:],         dtype=float)
        Mc    = np.asarray(d.variables["M_corrected"][:],    dtype=float)
        delta = np.asarray(d.variables["delta"][:],          dtype=float)
        cont  = np.asarray(d.variables["continent_mask"][:], dtype=bool)
    return lat, lon, M, Mc, delta, cont


# ---------------------------------------------------------------------------
# Sample overlay
# ---------------------------------------------------------------------------
_GEOCHEM_CACHE = None
def _get_geochem():
    """Cache the geochem CSV (with cached S&W plate IDs) across the three
    figures we render in one main() call."""
    global _GEOCHEM_CACHE
    if _GEOCHEM_CACHE is None:
        _GEOCHEM_CACHE = A.get_geochem()
    return _GEOCHEM_CACHE


def _declustered_samples_at(t, lat, lon, cont, M, Mc):
    """The per-slice declustered samples that actually entered the
    assimilation, reconstructed to age t.  Applies the strict
    visual-overlay land filter `cont & (M>0) & (Mc>0)` so points don't
    plot over rendered ocean."""
    df = A.prepare_samples(_get_geochem(), t)
    if df.empty:
        return df.assign(rlat=[], rlon=[], z_obs_m=[])
    iy = nearest_cell_index(lat, df["rlat"].values)
    ix = nearest_cell_index(lon, df["rlon"].values)
    on_land = cont[iy, ix] & (M[iy, ix] > 0) & (Mc[iy, ix] > 0)
    df = df[on_land].reset_index(drop=True)
    if df.empty:
        return df.assign(rlat=[], rlon=[], z_obs_m=[])
    dec = A.decluster(df, lat, lon)
    # decluster() emits columns: rlat, rlon, z (the weighted-median elevation), …
    return dec.rename(columns={"z": "z_obs_m"})


def panel(ax, lat, lon, arr, cmap, norm, cont, title, age_label, kind="elev",
          samples=None, coast_field=None):
    """Render one elevation/Δz panel.

    If `samples` is a DataFrame with rlon/rlat/z_obs_m columns it is overlaid:
      - on elevation panels: filled circles coloured by z_obs on the same
        terrain cmap as the background, so visual agreement is obvious.
      - on Δz panels: small open black circles marking sample positions only.

    `coast_field` is the elevation grid whose ≥0 contour is drawn as
    the coastline outline.  Pass M_orig (or M_corrected for the
    corrected-elevation panels) — this ensures the coastline matches
    the actual rendered topography rather than a separately-derived
    continental-polygon raster.
    """
    ax.set_global()
    LON2D, LAT2D = np.meshgrid(lon, lat)
    if kind == "delta":
        # Mask Δz to the actual paleo-land footprint of the rendered
        # field (`coast_field ≥ 0`) rather than the continent-polygon
        # raster, so we don't paint Δz on cells that sit inside the
        # polygon footprint but below sea level in S&W18.
        if coast_field is not None:
            land = (coast_field >= 0)
            arr = np.where(land, arr, np.nan)
        else:
            arr = np.where(cont, arr, np.nan)
    pcm = ax.pcolormesh(LON2D, LAT2D, arr, cmap=cmap, norm=norm,
                        transform=ccrs.PlateCarree(), shading="auto", rasterized=True)
    # Coastline outline: drawn ONLY on delta-z panels.  On elevation
    # panels the terrain colour ramp already makes the coastline visually
    # obvious, so the extra black contour is redundant clutter; on delta
    # panels the diverging colour scheme on a mostly-near-zero field
    # benefits from an explicit continental outline.
    if kind == "delta":
        coast = coast_field if coast_field is not None else arr
        ax.contour(LON2D, LAT2D, (coast >= 0).astype(float),
                   levels=[0.5], colors="black", linewidths=0.4,
                   transform=ccrs.PlateCarree())

    if samples is not None and len(samples):
        if kind == "elev":
            # Same colour mapping as the underlying terrain so the reader
            # can compare sample-derived z against the corrected field.
            ax.scatter(samples["rlon"].values, samples["rlat"].values,
                       c=samples["z_obs_m"].values,
                       cmap=plt.cm.terrain, norm=ELEV_NORM,
                       s=14, edgecolor="black", linewidths=0.35,
                       transform=ccrs.PlateCarree(), zorder=5)
        else:
            # Δz panel: position-only marker, no fill, doesn't fight the
            # diverging RdBu_r colour scheme.
            ax.scatter(samples["rlon"].values, samples["rlat"].values,
                       facecolors="none", edgecolor="black",
                       s=12, linewidths=0.5,
                       transform=ccrs.PlateCarree(), zorder=5)
        ax.text(0.02, 0.04, f"n = {len(samples)}",
                transform=ax.transAxes, fontsize=8,
                ha="left", va="bottom",
                bbox=dict(facecolor="white", alpha=0.7, edgecolor="none", pad=1.5))

    ax.gridlines(lw=0.25, color="grey", alpha=0.4)
    ax.set_title(title, fontsize=11)
    ax.text(0.97, 0.04, f"{int(age_label)} Ma",
            transform=ax.transAxes, fontsize=10, weight="bold",
            ha="right", va="bottom",
            bbox=dict(facecolor="white", alpha=0.8, edgecolor="none", pad=2))
    return pcm


def make_one(age_L, age_R, basename):
    print(f"  → {basename}: {age_L} Ma | {age_R} Ma")
    lat, lon, M_L, Mc_L, dL, cont_L = load(age_L)
    _, _,     M_R, Mc_R, dR, cont_R = load(age_R)

    # Per-slice declustered samples that entered the assimilation.
    samples_L = _declustered_samples_at(age_L, lat, lon, cont_L, M_L, Mc_L)
    samples_R = _declustered_samples_at(age_R, lat, lon, cont_R, M_R, Mc_R)
    print(f"     n_declustered: {age_L} Ma → {len(samples_L)},  "
          f"{age_R} Ma → {len(samples_R)}")

    fig = plt.figure(figsize=(13, 11.5))
    proj = ccrs.Robinson()
    pcm_elev = pcm_delta = None

    layout = [
        # (row, kind, arr_L, arr_R, coast_L, coast_R, title, show_samples)
        (0, "elev",  M_L,  M_R,  M_L,  M_R,  "Scotese & Wright 2018 input",       False),
        (1, "elev",  Mc_L, Mc_R, Mc_L, Mc_R, "Corrected elevation + samples",     True),
        (2, "delta", dL,   dR,   Mc_L, Mc_R, "Δz = corrected − original + samples", True),
    ]
    for row, kind, aL, aR, coastL, coastR, title_prefix, show_samples in layout:
        cmap = plt.cm.terrain if kind == "elev" else plt.cm.RdBu_r
        norm = ELEV_NORM if kind == "elev" else DELTA_NORM
        sL = samples_L if show_samples else None
        sR = samples_R if show_samples else None
        ax_L = fig.add_subplot(3, 2, row*2 + 1, projection=proj)
        ax_R = fig.add_subplot(3, 2, row*2 + 2, projection=proj)
        p_L = panel(ax_L, lat, lon, aL, cmap, norm, cont_L,
                    f"{title_prefix}", age_L, kind,
                    samples=sL, coast_field=coastL)
        p_R = panel(ax_R, lat, lon, aR, cmap, norm, cont_R,
                    f"{title_prefix}", age_R, kind,
                    samples=sR, coast_field=coastR)
        if kind == "elev":
            pcm_elev = p_L
        else:
            pcm_delta = p_L

    # Two colorbars at the bottom, one for elevation rows (1+2) and one for Δz row (3).
    cax_elev = fig.add_axes([0.18, 0.40, 0.64, 0.012])
    cb_elev = fig.colorbar(pcm_elev, cax=cax_elev, orientation="horizontal",
                            extend="both")
    cb_elev.set_label("elevation (m)", fontsize=10)
    cax_delta = fig.add_axes([0.18, 0.06, 0.64, 0.012])
    cb_delta = fig.colorbar(pcm_delta, cax=cax_delta, orientation="horizontal",
                             extend="both")
    cb_delta.set_label("Δz (m)", fontsize=10)

    fig.subplots_adjust(left=0.02, right=0.98, top=0.97, bottom=0.12,
                        hspace=0.35, wspace=0.05)
    out = FIG_DIR / f"{basename}_SW_comparison_{age_L}-{age_R}Ma.png"
    plt.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"     wrote {out}")


def main():
    for age_L, age_R, basename in PAIRS:
        make_one(age_L, age_R, basename)
    print("\nDone.  Three figures in paper/Scotese/.")


if __name__ == "__main__":
    main()
