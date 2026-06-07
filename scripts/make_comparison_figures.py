#!/usr/bin/env python3
"""
=============================================================================
make_comparison_figures.py  —  Multi-panel pyGMT figures comparing the
original Scotese & Wright PaleoDEMs against the geochemically-corrected
maps, with the Δz difference, side by side at selected ages.
=============================================================================

WHAT THIS DOES
    For each requested pair of ages, produces a 2-column × 3-row pyGMT
    figure in Winkel-Tripel projection.  Columns = ages, rows =
        row 1 — original Scotese & Wright elevation
        row 2 — geochemically-corrected elevation
        row 3 — Δz = corrected − original.
    Shared colorbars at the bottom (one for elevation rows, one for Δz).

    Default age pairs (one figure per pair):
        figure 1: (500 Ma, 400 Ma)
        figure 2: (300 Ma, 200 Ma)
        figure 3: (100 Ma,  50 Ma)

OUTPUT
    paths_scotese.OUTPUT_DIR / "SW_comparison_<age1>-<age2>Ma.png"
    paths_scotese.OUTPUT_DIR / "SW_comparison_<age1>-<age2>Ma.pdf"

USAGE
    cd <project>/scripts
    python make_comparison_figures.py                       # default 6 ages
    python make_comparison_figures.py --pairs 50 100 200 400 # 50&100, 200&400
    python make_comparison_figures.py --no-pdf              # PNG only
    python make_comparison_figures.py --dpi 300             # publication PNG

OPTIONS
    --pairs  AGE1 AGE2 [AGE3 AGE4 ...]
        Explicit age list — even number of integers; consecutive pairs become
        figures.  Default: 500 400 300 200 100 0.
    --no-pdf            skip PDF export
    --no-png            skip PNG export
    --dpi  INT          PNG dpi (default 200)
    --width-cm FLOAT    width of each map panel (default 10 cm)

DEPENDENCIES
    GMT 6.x, pygmt, xarray, netCDF4
    Optional: pygplates (only for subduction-zone overlay where resolvable —
                        0–100 Ma in S&W).  Set DRAW_SZ=False to skip.

NOTES
    - cpt range for elevation is FIXED at ±4000 m with end-of-bar triangle
      arrows (`+e`) so the same colour scale is used in every panel of every
      figure for direct visual comparison.
    - cpt range for Δz is FIXED at ±2000 m (= the assimilation's per-cell cap).
    - This script intentionally calls pyGMT for **all** plotting so that the
      output is publication-grade Winkel-Tripel, matching the videos.
=============================================================================
"""
from __future__ import annotations
import argparse, sys, shutil
from pathlib import Path
import numpy as np
import netCDF4 as nc
import xarray as xr
import pygmt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from paths_scotese import CORRECTED_DIR, OUTPUT_DIR, PROJECT_ROOT
from sw_io import available_ages as sw_available_ages, nearest_cell_index
import assimilate_scotese as A   # for prepare_samples + decluster


# ---------------------------------------------------------------------------
# Geochem sample overlay (mirrors the cartopy Fig 5 + delta-mode videos)
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
    """Per-slice declustered samples that actually entered the
    assimilation, reconstructed to age `t`.  Applies the strict
    visual-overlay land filter `cont & (M>0) & (Mc>0)` so points don't
    plot over rendered ocean.  Returns an empty DataFrame if none.
    """
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
    return dec.rename(columns={"z": "z_obs_m"})

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
ELEV_CPT = "earth"
DELTA_CPT = "polar"
ELEV_RANGE = (-4000.0, 4000.0, 250.0)     # series for makecpt
DELTA_RANGE = (-2000.0, 2000.0, 100.0)
PROJ_BASE = "R0"                          # Winkel-Tripel centred at lon=0
REGION = "g"                              # global
DRAW_SZ = True                            # overlay subduction zones where resolvable

# Default age pairs (consecutive pairs of the list become figures)
DEFAULT_AGES = [500, 400, 300, 200, 100, 50]


# ---------------------------------------------------------------------------
# Data helpers
# ---------------------------------------------------------------------------
def load_grid(t: int):
    """Return (M_orig, M_corrected, delta, lat, lon, continent_mask)."""
    f = CORRECTED_DIR / f"{int(t)}Ma_corrected_SW.nc"
    if not f.exists():
        raise FileNotFoundError(
            f"{f} not found — run `python assimilate_scotese.py {t}` first")
    with nc.Dataset(f) as d:
        M     = d.variables["M_orig"][:].astype(float)
        Mc    = d.variables["M_corrected"][:].astype(float)
        delta = d.variables["delta"][:].astype(float)
        lat   = d.variables["lat"][:].astype(float)
        lon   = d.variables["lon"][:].astype(float)
        cont  = d.variables["continent_mask"][:].astype(bool)
    return M, Mc, delta, lat, lon, cont


def to_xr(arr: np.ndarray, lat, lon) -> xr.DataArray:
    """Wrap a (lat, lon) ndarray as an xarray DataArray, tagged as
    pixel-registered + geographic so GMT applies the correct periodic
    dateline wrap and doesn't emit the "Longitude range too small;
    geographic boundary condition changed to natural" warning.

    The S&W corrected NetCDFs are pixel-registered (cell-centre lons
    run -179.5..+179.5, so centre-to-centre span is 359° even though
    the underlying coverage is 360°).  Without these attrs GMT
    defaults to gridline-registered Cartesian, reads the 359° span
    literally, and skips the periodic wrap.
    """
    da = xr.DataArray(arr.astype(np.float32),
                      coords={"lat": lat, "lon": lon},
                      dims=("lat", "lon"), name="z")
    da.gmt.registration = 1   # 0 = gridline, 1 = pixel
    da.gmt.gtype = 1          # 0 = cartesian, 1 = geographic
    return da


def subduction_zones(t):
    """Return list of (lat, lon) arrays — empty if topologies don't resolve at t."""
    if not DRAW_SZ:
        return []
    try:
        from plate_model_utils_scotese import subduction_zones as _sz
        return _sz(float(t))
    except Exception as e:
        print(f"  subduction zone overlay skipped at {t} Ma: {e}")
        return []


# ---------------------------------------------------------------------------
# Single panel drawer (one map within the multi-panel figure)
# ---------------------------------------------------------------------------
def draw_panel(fig, ds_da: xr.DataArray, cpt_path: str, title: str,
               proj: str, region: str = REGION,
               cont_da: xr.DataArray = None,
               sz_lines = None,
               show_age_label: str = None,
               show_coast: bool = True,
               samples=None):
    """Render one map.  Caller is responsible for shifting origin.

    show_coast : draw the gray coastline contour derived from `cont_da`.
                 Switched off for the elevation rows where the terrain
                 colour ramp already makes the coastline visually
                 obvious; left on for the delta-z row where the diverging
                 colour scheme on a mostly-near-zero field benefits from
                 an explicit continental outline.

    samples    : optional DataFrame with rlon/rlat columns — the
                 declustered geochemical samples that drove the
                 assimilation at this age.  Drawn as small open black
                 circles (position-only marker; z-colouring on a
                 diverging Δz cmap would mislead because z_obs is
                 absolute and Δz is a signed correction).  Mirrors the
                 delta-mode overlay in the publication videos.
    """
    # Note: the title is passed as `+t<title>` WITHOUT literal quotes —
    # adding inner double quotes (`+t"..."`) causes them to render
    # visibly around the title in the output PNG/PDF.
    fig.basemap(region=region, projection=proj, frame=["af", f"+t{title}"])
    fig.grdimage(grid=ds_da, projection=proj, region=region,
                 cmap=cpt_path, nan_transparent=True)
    if show_coast and cont_da is not None and cont_da.values.any():
        fig.grdcontour(grid=cont_da, projection=proj, region=region,
                       levels=[0.5], pen="0.4p,gray30")
    if sz_lines:
        for line in sz_lines:
            xy = np.column_stack([line[:, 1], line[:, 0]])
            fig.plot(x=xy[:, 0], y=xy[:, 1],
                     pen="0.5p,red", projection=proj, region=region)
    if samples is not None and len(samples):
        fig.plot(x=samples["rlon"].values, y=samples["rlat"].values,
                 style="c0.10c", pen="0.5p,black",
                 projection=proj, region=region)
    if show_age_label:
        # Upper-left of the map.  Anchor is in lon/lat; the cm offset is
        # applied AFTER projection, so it's predictable regardless of
        # how Winkel-Tripel squashes the meridians at high latitude.
        # Tweak the offset (left/up in cm) to nudge the label.
        fig.text(x=-170, y=65, text=show_age_label,
                 offset="-1.8c/0.8c", # 1.8 cm left, 0.8 cm up from anchor
                 font="14p,Helvetica-Bold,black",
                 region=region, projection=proj, justify="ML",
                 no_clip=True)


# ---------------------------------------------------------------------------
# Build one 2 × 3 comparison figure
# ---------------------------------------------------------------------------
def make_figure(age_left: int, age_right: int, out_basename: str,
                width_cm: float = 10.0, write_png: bool = True,
                write_pdf: bool = True, dpi: int = 200,
                paper_basename: str | None = None):
    """
    paper_basename : if given, also write the figure to
        PROJECT_ROOT/paper/Scotese/<paper_basename>.{png,pdf}
        in addition to OUTPUT_DIR/<out_basename>.{png,pdf}.  Used to
        publish the canonical Fig05a/b/c renderings without an extra
        shell-side copy.
    """
    print(f"\n=== building comparison figure: {age_left} Ma | {age_right} Ma ===")
    proj = f"{PROJ_BASE}/{width_cm}c"  # e.g. "R0/10c"

    # Load both ages
    Ma_L, Mc_L, dL, lat, lon, cont_L = load_grid(age_left)
    Ma_R, Mc_R, dR, lat2, lon2, cont_R = load_grid(age_right)
    assert lat.shape == lat2.shape and lon.shape == lon2.shape, "grid mismatch"

    sz_L = subduction_zones(age_left)
    sz_R = subduction_zones(age_right)

    # Coastline contour rasters derived from the actual rendered
    # elevation field (M ≥ 0).  Two flavours: one for the "original"
    # row (rendered M_orig), one for the corrected + Δz rows (rendered
    # M_corrected).  This avoids the polygon vs S&W18 paleo-coastline
    # mismatch that becomes visible before ~350 Ma.
    coast_L_orig_da = to_xr((Ma_L  >= 0).astype(np.float32), lat, lon)
    coast_R_orig_da = to_xr((Ma_R  >= 0).astype(np.float32), lat, lon)
    coast_L_corr_da = to_xr((Mc_L  >= 0).astype(np.float32), lat, lon)
    coast_R_corr_da = to_xr((Mc_R  >= 0).astype(np.float32), lat, lon)

    # Per-slice declustered samples that actually entered the
    # assimilation, reconstructed to each target age.  Overlaid on the
    # delta-z row as position markers (same convention as the videos +
    # the cartopy Fig 5 fallback).
    samples_L = _declustered_samples_at(age_left,  lat, lon, cont_L, Ma_L, Mc_L)
    samples_R = _declustered_samples_at(age_right, lat, lon, cont_R, Ma_R, Mc_R)
    print(f"     n_declustered: {age_left} Ma -> {len(samples_L)},  "
          f"{age_right} Ma -> {len(samples_R)}")

    # Build the two CPT files
    elev_cpt = str(OUTPUT_DIR / "_tmp_elev.cpt")
    delta_cpt = str(OUTPUT_DIR / "_tmp_delta.cpt")
    pygmt.makecpt(cmap=ELEV_CPT, series=ELEV_RANGE,
                  continuous=True, background=True, output=elev_cpt)
    pygmt.makecpt(cmap=DELTA_CPT, series=DELTA_RANGE,
                  continuous=True, background=True, output=delta_cpt)

    fig = pygmt.Figure()
    pygmt.config(MAP_FRAME_TYPE="plain",
                 FONT_TITLE="12p,Helvetica-Bold,black",
                 FONT_ANNOT_PRIMARY="8p,Helvetica,black",
                 FONT_LABEL="9p,Helvetica,black",
                 MAP_TITLE_OFFSET="-2p")   # ~1.4 mm tighter than default 2p

    # Approximate panel height for Winkel-Tripel at this width:
    # the Winkel-Tripel aspect ratio is roughly 1.637:1 (width:height).
    panel_h = width_cm / 1.637
    h_gap = 1.4        # vertical gap between rows (in cm)
    v_gap = 0.6        # horizontal gap between cols (in cm)
    # Extra space between row 2 (corrected) and row 3 (delta z) to host
    # the shared elevation colorbar — without this, the colorbar overlays
    # the bottom of the Corrected maps.  Tuned so the colorbar sits ~1.4 cm
    # below row 2 (see the colorbar's `-1.4c` y-offset further down) and
    # row 3 is comfortably below the colorbar.
    cb_band = 1.4      # cm; height of the elevation-colorbar band

    # Origin starts at bottom-left; we draw bottom row first to top row.
    # That way subsequent shifts ADD to y as we go up.

    # ---- ROW 3 (bottom): delta z ----
    # ASCII-only titles + colorbar label — Ghostscript 10.x doesn't
    # define the ISOLatin1+_Encoding GMT emits when non-ASCII glyphs
    # appear in PostScript, so we spell out Greek delta as "delta",
    # true minus as "-", plus-minus as "+/-".  Matplotlib-rendered
    # figures elsewhere keep the Unicode versions.
    fig.shift_origin(xshift="0c", yshift="3.5c")    # leave room for colorbars
    # left panel
    draw_panel(fig, to_xr(dL, lat, lon), delta_cpt,
               title="delta z",
               proj=proj, cont_da=coast_L_corr_da, sz_lines=sz_L,
               show_age_label=f"{age_left} Ma",
               samples=samples_L)
    # right panel
    fig.shift_origin(xshift=f"{width_cm+v_gap}c")
    draw_panel(fig, to_xr(dR, lat, lon), delta_cpt,
               title="delta z",
               proj=proj, cont_da=coast_R_corr_da, sz_lines=sz_R,
               show_age_label=f"{age_right} Ma",
               samples=samples_R)
    # delta z colorbar across both bottom panels
    cb_total_w = 2*width_cm + v_gap
    fig.colorbar(projection=proj, region=REGION, cmap=delta_cpt,
                 frame=['x+l"delta z  corrected - original  (m)"', "af"],
                 position=f"jBC+w{cb_total_w*0.6:.1f}c/0.3c+o-{(width_cm+v_gap)/2:.1f}c/-1.6c+h+e")

    # back to the left column for the middle row.  The y-shift is bigger
    # than the normal inter-row gap so we have a `cb_band`-tall slot
    # between rows 2 and 3 that hosts the shared elevation colorbar.
    fig.shift_origin(xshift=f"-{width_cm+v_gap}c",
                     yshift=f"{panel_h+h_gap+cb_band}c")

    # ---- ROW 2 (middle): corrected ----
    draw_panel(fig, to_xr(Mc_L, lat, lon), elev_cpt,
               title="Corrected",
               proj=proj, cont_da=coast_L_corr_da, sz_lines=sz_L,
               show_age_label=f"{age_left} Ma",
               show_coast=False)
    fig.shift_origin(xshift=f"{width_cm+v_gap}c")
    draw_panel(fig, to_xr(Mc_R, lat, lon), elev_cpt,
               title="Corrected",
               proj=proj, cont_da=coast_R_corr_da, sz_lines=sz_R,
               show_age_label=f"{age_right} Ma",
               show_coast=False)

    # elevation colorbar — placed immediately under row 2, centered across
    # both columns, sitting in the `cb_band` gap we left between rows 2
    # and 3.  Mirrors the delta-z colorbar layout: x-offset is half a
    # column-pitch to the LEFT of `jBC` (anchored on the right panel),
    # which lands it under the boundary between the two panels; y-offset
    # is downward into the gap.
    fig.colorbar(projection=proj, region=REGION, cmap=elev_cpt,
                 frame=['x+l"elevation (m)"', "af"],
                 position=f"jBC+w{cb_total_w*0.6:.1f}c/0.3c+o-{(width_cm+v_gap)/2:.1f}c/-1.4c+h+e")

    # ---- ROW 1 (top): original ----
    fig.shift_origin(xshift=f"-{width_cm+v_gap}c",
                     yshift=f"{panel_h+h_gap}c")
    draw_panel(fig, to_xr(Ma_L, lat, lon), elev_cpt,
               title="Scotese and Wright (2018)",
               proj=proj, cont_da=coast_L_orig_da, sz_lines=sz_L,
               show_age_label=f"{age_left} Ma",
               show_coast=False)
    fig.shift_origin(xshift=f"{width_cm+v_gap}c")
    draw_panel(fig, to_xr(Ma_R, lat, lon), elev_cpt,
               title="Scotese and Wright (2018)",
               proj=proj, cont_da=coast_R_orig_da, sz_lines=sz_R,
               show_age_label=f"{age_right} Ma",
               show_coast=False)

    # ---- Save ----
    # Each output is written to two locations:
    #   1) OUTPUT_DIR (outputs_Scotese/SW_comparison_*.{png,pdf}) — the
    #      intermediate render path; cheap to wipe and rebuild during
    #      layout iteration.
    #   2) paper/Scotese/ with the publication name (only when
    #      `paper_basename` is given) — so a direct `python
    #      make_comparison_figures.py` run refreshes the paper-folder
    #      copy in one step, no extra `cp` needed.
    paper_dir = PROJECT_ROOT / "paper" / "Scotese"
    if paper_basename:
        paper_dir.mkdir(parents=True, exist_ok=True)

    if write_png:
        out_png = OUTPUT_DIR / f"{out_basename}.png"
        if out_png.exists(): out_png.unlink()
        fig.savefig(out_png, dpi=dpi)
        print(f"  wrote {out_png}")
        if paper_basename:
            paper_png = paper_dir / f"{paper_basename}.png"
            if paper_png.exists(): paper_png.unlink()
            fig.savefig(paper_png, dpi=dpi)
            print(f"  wrote {paper_png}")
    if write_pdf:
        out_pdf = OUTPUT_DIR / f"{out_basename}.pdf"
        if out_pdf.exists(): out_pdf.unlink()
        fig.savefig(out_pdf)
        print(f"  wrote {out_pdf}")
        if paper_basename:
            paper_pdf = paper_dir / f"{paper_basename}.pdf"
            if paper_pdf.exists(): paper_pdf.unlink()
            fig.savefig(paper_pdf)
            print(f"  wrote {paper_pdf}")

    # Clean up temporary CPTs
    for f in [elev_cpt, delta_cpt]:
        try: Path(f).unlink()
        except FileNotFoundError: pass


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pairs", type=int, nargs="+", default=DEFAULT_AGES,
                   help="ages to compare, in pairs (e.g. 500 400 300 200 100 0)")
    p.add_argument("--no-png", action="store_true")
    p.add_argument("--no-pdf", action="store_true")
    p.add_argument("--dpi", type=int, default=200)
    p.add_argument("--width-cm", type=float, default=10.0,
                   help="width of each map panel in cm (default 10)")
    args = p.parse_args()

    ages = args.pairs
    if len(ages) % 2 != 0:
        raise SystemExit("ERROR: --pairs must contain an even number of ages")

    # Validate that requested ages have corrected NetCDFs
    have = set(sw_available_ages())
    missing = [a for a in ages if a not in have]
    if missing:
        raise SystemExit(f"ERROR: no S&W slice at ages {missing}.  Available: {sorted(have)}")

    # Canonical paper-figure mapping for the three age pairs used in
    # Earth-Science Reviews Fig 5.  Any other pair is rendered only to
    # OUTPUT_DIR (no paper-folder copy).
    PAPER_PAIRS = {
        (500, 400): "Fig05a_SW_comparison_500-400Ma",
        (300, 200): "Fig05b_SW_comparison_300-200Ma",
        (100,  50): "Fig05c_SW_comparison_100-50Ma",
    }

    for i in range(0, len(ages), 2):
        a1, a2 = ages[i], ages[i+1]
        out_name = f"SW_comparison_{a1}-{a2}Ma"
        paper_name = PAPER_PAIRS.get((a1, a2))   # None for non-canonical pairs
        make_figure(a1, a2, out_basename=out_name,
                    paper_basename=paper_name,
                    width_cm=args.width_cm,
                    write_png=not args.no_png,
                    write_pdf=not args.no_pdf,
                    dpi=args.dpi)


if __name__ == "__main__":
    main()
