"""
=============================================================================
make_pipeline_illustration_figure.py  —  Figure 3 of the Earth-Science Reviews paper
=============================================================================

Six-panel illustration of the assimilation pipeline at one representative
time slice (100 Ma):

  (a) Scotese & Wright 2018 input elevation
  (b) Reconstructed declustered samples (filled circles, coloured by
      observed elevation) overlaid on the input grid
  (c) Tectonic-province raster with subduction-zone overlay
  (d) Δz field from province-wise CDF rescaling only (after the 400 km
      support gate, the province-edge taper and the cross-boundary smoothing)
  (e) Δz field from smoothed-residual kernel only (after its 250-500 km
      reach taper)
  (f) Final corrected map (after continent-masked smoothing + cap)

Every panel is drawn from the intermediate fields of ONE call to
assimilate_scotese.assimilate_one (its `stages` argument), so the figure
shows the chain the product is built with; panel (f) is checked against
the product grid of the same age in CORRECTED_DIR before anything is drawn.

All panels in Robinson projection.  Elevation cpt: terrain, −5..+5 km.
Δz cpt: RdBu_r, ±2 km.  Province cmap: tab10.

OUTPUT
    Figures/Fig03_pipeline_illustration_100Ma.png

USAGE
    cd <project>/scripts
    python make_pipeline_illustration_figure.py
=============================================================================
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import netCDF4 as nc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from paths_scotese import CORRECTED_DIR, PAPER_FIG_DIR
import assimilate_scotese as A
from plate_model_utils_scotese import subduction_zones

PROJ_ROOT = HERE.parent
# Paper-numbered output: land directly in Figures/, matching the
# convention used by the other Fig0X scripts in this directory.
FIG_DIR = PAPER_FIG_DIR
FIG_DIR.mkdir(parents=True, exist_ok=True)

AGE = 100
ELEV_NORM = mcolors.TwoSlopeNorm(vmin=-5000, vcenter=0, vmax=5000)
DELTA_NORM = mcolors.TwoSlopeNorm(vmin=-2000, vcenter=0, vmax=2000)


def add_panel(fig, idx, projection):
    return fig.add_subplot(3, 2, idx, projection=projection)


def show_topo(ax, lat, lon, arr, cont, title, label_letter, label_age=True,
              cmap=plt.cm.terrain, norm=ELEV_NORM, kind="elev",
              coast_field=None):
    """`coast_field` is the elevation grid whose ≥0 contour is drawn.
    Defaults to `arr` for elevation panels; on Δz panels (d/e) pass the
    underlying M so the contour traces the actual coastline rather
    than the Δz=0 isoline."""
    ax.set_global()
    LON2D, LAT2D = np.meshgrid(lon, lat)
    if kind == "delta":
        arr = np.where(cont, arr, np.nan)
    pcm = ax.pcolormesh(LON2D, LAT2D, arr, cmap=cmap, norm=norm,
                        transform=ccrs.PlateCarree(), shading="auto", rasterized=True)
    coast = coast_field if coast_field is not None else arr
    ax.contour(LON2D, LAT2D, (coast >= 0).astype(float),
               levels=[0.5], colors="black", linewidths=0.4,
               transform=ccrs.PlateCarree())
    ax.gridlines(lw=0.25, color="grey", alpha=0.4)
    ax.set_title(f"{label_letter}  {title}", fontsize=11, loc="left")
    if label_age:
        ax.text(0.97, 0.05, f"{int(AGE)} Ma", transform=ax.transAxes,
                fontsize=9, weight="bold", ha="right", va="bottom",
                bbox=dict(facecolor="white", alpha=0.8, edgecolor="none", pad=2))
    return pcm


def main():
    print(f"Rendering pipeline illustration at {AGE} Ma …")
    # One run of the workflow for this slice, keeping its intermediate fields, so each
    # panel is exactly the stage the product went through (the figure used to recompute the
    # stages itself and left out the support gate, the edge taper, the cross-boundary
    # smoothing and the reach taper).
    st = {}
    A.assimilate_one(float(AGE), save_nc=False, stages=st)
    M, cont, dec, P, M3, r, Mf = (st[k] for k in ("M", "cont", "dec", "P", "M3", "r", "M_final"))
    _, lat, lon = A.sw_load_grid(AGE)
    prod = CORRECTED_DIR / f"{AGE}Ma_corrected_SW.nc"
    with nc.Dataset(prod) as d:
        Mp = np.asarray(d["M_corrected"][:], dtype=float)
    dmax = float(np.nanmax(np.abs(Mf - Mp)))
    print(f"  panel (f) against the product {prod.name}: max |difference| {dmax:.4f} m")
    if dmax > 0.5:          # exact on the machine that built the product; other numpy/scipy builds differ by < 0.1 m
        sys.exit(f"  the workflow run does not reproduce {prod}: the figure is not drawn")

    fig = plt.figure(figsize=(14, 13))
    proj = ccrs.Robinson()

    # (a) input
    ax = add_panel(fig, 1, proj)
    pcm_elev = show_topo(ax, lat, lon, M, cont,
                         "Scotese & Wright 2018 input", "(a)")

    # (b) samples on top of the input map
    ax = add_panel(fig, 2, proj)
    show_topo(ax, lat, lon, M, cont,
              f"Declustered samples (n={len(dec)})", "(b)")
    ax.scatter(dec["rlon"], dec["rlat"], c=dec["z"], cmap="terrain", norm=ELEV_NORM,
               s=22, edgecolor="black", linewidths=0.5,
               transform=ccrs.PlateCarree(), zorder=5)

    # (c) province grid with SZ overlay
    ax = add_panel(fig, 3, proj)
    ax.set_global()
    LON2D, LAT2D = np.meshgrid(lon, lat)
    P_plot = np.where(cont, P.astype(float), np.nan)
    ax.pcolormesh(LON2D, LAT2D, P_plot, cmap="tab10", vmin=0, vmax=9,
                  transform=ccrs.PlateCarree(), shading="auto", rasterized=True)
    ax.contour(LON2D, LAT2D, (M >= 0).astype(float),
               levels=[0.5], colors="black", linewidths=0.4,
               transform=ccrs.PlateCarree())
    ax.gridlines(lw=0.25, color="grey", alpha=0.4)
    # Subduction zones overlay
    for line in subduction_zones(float(AGE)):
        ax.plot(line[:, 1], line[:, 0], color="red", lw=0.8,
                transform=ccrs.PlateCarree())
    ax.set_title("(c)  Tectonic-province raster (SZs in red)", fontsize=11, loc="left")
    ax.text(0.97, 0.05, f"{AGE} Ma", transform=ax.transAxes,
            fontsize=9, weight="bold", ha="right", va="bottom",
            bbox=dict(facecolor="white", alpha=0.8, edgecolor="none", pad=2))

    # (d) Δz from CDF rescaling only
    ax = add_panel(fig, 4, proj)
    pcm_d = show_topo(ax, lat, lon, M3 - M, cont,
                      "Δz from province CDF rescaling only", "(d)",
                      cmap=plt.cm.RdBu_r, norm=DELTA_NORM, kind="delta",
                      coast_field=M)

    # (e) Δz from residual kernel only
    ax = add_panel(fig, 5, proj)
    show_topo(ax, lat, lon, r, cont,
              "Δz from smoothed-residual kernel only", "(e)",
              cmap=plt.cm.RdBu_r, norm=DELTA_NORM, kind="delta",
              coast_field=M)

    # (f) final corrected
    ax = add_panel(fig, 6, proj)
    show_topo(ax, lat, lon, Mf, cont,
              "Final corrected (smoothed + capped)", "(f)")

    # Shared colorbars
    cax_e = fig.add_axes([0.05, 0.04, 0.42, 0.012])
    cb_e = fig.colorbar(pcm_elev, cax=cax_e, orientation="horizontal", extend="both")
    cb_e.set_label("elevation (m)", fontsize=9)
    cax_d = fig.add_axes([0.55, 0.04, 0.42, 0.012])
    cb_d = fig.colorbar(pcm_d, cax=cax_d, orientation="horizontal", extend="both")
    cb_d.set_label("Δz (m)", fontsize=9)

    fig.subplots_adjust(left=0.02, right=0.98, top=0.94, bottom=0.10,
                         hspace=0.28, wspace=0.06)
    out = FIG_DIR / "Fig03_pipeline_illustration_100Ma.png"
    plt.savefig(out, dpi=180, bbox_inches="tight")
    print(f"  wrote {out}")


if __name__ == "__main__":
    main()
