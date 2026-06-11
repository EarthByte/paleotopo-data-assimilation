#!/usr/bin/env python3
"""
=============================================================================
make_dyntopo_panels_figure.py  —  Supp-mat 6-panel overview of the Young
                                  2022 dyntopo PER-STEP INCREMENT (plate-frame
                                  computed, rotated to Scotese paleomag frame)
=============================================================================

WHAT THIS DOES
    Produces a single 2-row x 3-col pyGMT figure in Winkel-Tripel
    projection showing the Young 2022 dyntopo PER-STEP INCREMENT
    z_dyntopo_diff = [dyntopo(t) - dyntopo(t - Δt)],  Δt = 5 Myr by default at six fixed ages.
    Each panel shows the same correction signal that
    `build_dyntopo_diff_correction.py` adds to M_corrected at that age,
    laid out as a stand-alone 2x3 supp-mat figure.

    Default layout: 300 / 250 / 200 Ma on the top row, 150 / 100 / 50 Ma
    on the bottom row.  The deepest age is 300 Ma rather than the 500
    Ma starting point of the main-text Fig 5 because the Young 2022
    GLD428 model output only extends to 300 Ma.  Polar diverging
    colourmap, single shared colourbar below the bottom row.

WHY THE PER-STEP INCREMENT (and not the cumulative anomaly)
    The present-day observed topography already contains today's
    dynamic-topography contribution.  Composing absolute past dyntopo
    onto a present-day-anchored reconstruction would double-count
    today's contribution.  The actual correction signal applied in
    `build_dyntopo_diff_correction.py` is the difference
        Δz(t) = dyntopo(t) - dyntopo(t - Δt),  Δt = 5 Myr by default
    which is zero at t=0 by construction (no predecessor).  The subtraction is performed
    in PLATE reference frame (each grid cell rigidly attached to its
    continent across time) so the comparison is between the same point
    on the same continent at two different times; subtracting at the
    same lat/lon in any non-plate-fixed frame would compare different
    parts of different continents and is physically meaningless.  The
    plate-frame difference is then cookie-cut by Scotese 2023
    continental polygons and rotated into the time-t Scotese paleomag
    frame for visualization — which is what's read here from
    `<age>Ma_corrected_plus_dyntopo_diff_young_SW.nc / z_dyntopo_diff`.

    This figure is intended for the supplementary material of the
    geochem-assimilation paper.  It illustrates the spatial and
    temporal pattern of the Young 2022 dyntopo field that informs the
    related supp-mat figure where dyntopo is additively composed with
    the geochem-corrected S&W maps (`make_comparison_figures_dyntopo.py`).

All six panel ages must have a combined NetCDF in --combined-dir
    (the output of build_dyntopo_diff_correction.py).  The script
    refuses to run if any are missing.

OUTPUT
    paths_scotese.OUTPUT_DIR / "dyntopo_diff_panels_300-50Ma.png"
    paths_scotese.OUTPUT_DIR / "dyntopo_diff_panels_300-50Ma.pdf"
    (filename is suffixed with the actual --ages list at run time)

USAGE
    cd <project>/scripts_Scotese
    # 1. (one-off) build the diff NetCDFs at the six target ages
    python build_dyntopo_diff_correction.py --source young \\
        --ages 50 100 150 200 250 300
    # 2. render the supp-mat panel figure
    python make_dyntopo_panels_figure.py

    # Custom age list (must be 6 ages, will be laid out top-row-first)
    python make_dyntopo_panels_figure.py --ages 250 200 150 100 50 0

OPTIONS
    --combined-dir PATH              directory of
                                     <age>Ma_corrected_plus_dyntopo_diff_young_SW.nc
                                     files (default:
                                     PROJECT_ROOT/data/corrected_Scotese_plus_dyntopo_diff_young)
    --ages INT INT INT INT INT INT   exactly six ages, top-row-first
                                     (default: 300 250 200 150 100 50)
    --vlim FLOAT                     symmetric +/- m for the polar cpt
                                     (default 500)
    --width-cm FLOAT                 width of each panel in cm (default 7)
    --no-pdf                         skip PDF export
    --no-png                         skip PNG export
    --dpi INT                        PNG dpi (default 200)

DEPENDENCIES
    GMT 6.x, pygmt, xarray, netCDF4

NOTES
    - cpt is GMT's `polar` (blue / white / red), symmetric +/- `--vlim` m.
    - Coastlines / continent outlines are NOT drawn: in the mantle ref
      frame the standard paleo-coastline polygons would need a separate
      mantle-frame rotation to align with the dyntopo.  The figure is
      intentionally the field-only view.
    - Each panel carries an in-map age label in the upper-left corner,
      matching the convention of make_comparison_figures.py.
=============================================================================
"""
from __future__ import annotations
import argparse, sys
from pathlib import Path
import numpy as np
import netCDF4 as nc
import xarray as xr
import pygmt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from paths_scotese import OUTPUT_DIR


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
DYNTOPO_CPT = "polar"
DEFAULT_VLIM = 500.0
PROJ_BASE = "R0"                  # Winkel-Tripel centred at lon=0
REGION = "g"                      # global
DEFAULT_AGES = [300, 250, 200, 150, 100, 50]   # top-row-first
# Note: Young 2022 GLD428 only extends to 300 Ma, so the deepest panel
# is 300 Ma rather than the 500 Ma that the main text figure starts at.
COMBINED_FNAME_FMT  = "{age}Ma_corrected_plus_dyntopo_diff_young_SW.nc"
DEFAULT_COMBINED_DIR = Path(__file__).resolve().parent.parent / "data" \
    / "corrected_Scotese_plus_dyntopo_diff_young"


# ---------------------------------------------------------------------------
# Data helpers
# ---------------------------------------------------------------------------
def load_dyntopo_diff(combined_dir: Path, age: int):
    """Return (z_dyntopo_diff, lat, lon) for the Scotese-paleomag-frame
    dyntopo per-step increment at `age`, read from the build_dyntopo_diff_
    correction.py output NetCDF.
    """
    f = combined_dir / COMBINED_FNAME_FMT.format(age=int(age))
    if not f.exists():
        raise FileNotFoundError(
            f"{f} not found - run "
            f"`python build_dyntopo_diff_correction.py --source young "
            f"--ages {int(age)}` first")
    with nc.Dataset(f) as d:
        lat = np.asarray(d["lat"][:], dtype=np.float64)
        lon = np.asarray(d["lon"][:], dtype=np.float64)
        z   = np.asarray(d["z_dyntopo_diff"][:], dtype=np.float32)
    return z, lat, lon


def to_xr(arr: np.ndarray, lat, lon) -> xr.DataArray:
    """Wrap a (lat, lon) ndarray as an xarray DataArray.

    Tagged as pixel-registered + geographic by default - this matches
    the convention used elsewhere in the project and silences GMT's
    "Longitude range too small; geographic boundary condition changed
    to natural" warning.  If your dyntopo grids turn out to be
    gridline-registered (lon = -180 .. +180 inclusive instead of
    -179.75 .. +179.75 at 0.5 deg), set registration=0 here.
    """
    da = xr.DataArray(arr.astype(np.float32),
                      coords={"lat": lat, "lon": lon},
                      dims=("lat", "lon"), name="z")
    da.gmt.registration = 1   # 0 = gridline, 1 = pixel
    da.gmt.gtype = 1          # 0 = cartesian, 1 = geographic
    return da


# ---------------------------------------------------------------------------
# Single panel drawer
# ---------------------------------------------------------------------------
def draw_panel(fig, ds_da: xr.DataArray, cpt_path: str,
               proj: str, region: str = REGION,
               age_label: str = None):
    """Render one panel.  Caller is responsible for shifting origin."""
    # Bare axis frame, no title at the top of each panel (the age
    # label sits inside the upper-left of the map instead).
    fig.basemap(region=region, projection=proj, frame="af")
    fig.grdimage(grid=ds_da, projection=proj, region=region,
                 cmap=cpt_path, nan_transparent=True)
    if age_label:
        # Upper-left of the map.  Anchor is in lon/lat; the cm offset
        # is applied AFTER projection, so it's predictable regardless
        # of how Winkel-Tripel squashes the meridians at high latitude.
        fig.text(x=-170, y=65, text=age_label,
                 offset="-1.8c/0.8c",
                 font="14p,Helvetica-Bold,black",
                 region=region, projection=proj, justify="ML",
                 no_clip=True)


# ---------------------------------------------------------------------------
# Build the 2 x 3 figure
# ---------------------------------------------------------------------------
def make_figure(ages: list[int], combined_dir: Path,
                out_basename: str,
                vlim: float = DEFAULT_VLIM,
                width_cm: float = 7.0, write_png: bool = True,
                write_pdf: bool = True, dpi: int = 200):
    assert len(ages) == 6, "exactly 6 ages required (3 cols x 2 rows)"
    print(f"\n=== dyntopo time-difference overview "
          f"(Scotese paleomag frame): "
          f"top row {ages[0:3]}, bottom row {ages[3:6]} ===")

    proj = f"{PROJ_BASE}/{width_cm}c"

    # Pull z_dyntopo_diff from each of the six combined NetCDFs.  The
    # plate-frame subtraction + cookie-cut + rotation has already
    # happened in build_dyntopo_diff_correction.py.
    panels = []
    for a in ages:
        z, lat, lon = load_dyntopo_diff(combined_dir, a)
        panels.append((a, z, lat, lon))
        print(f"  loaded {a} Ma: shape={z.shape}, "
              f"z_dyntopo_diff range=[{np.nanmin(z):.0f}, "
              f"{np.nanmax(z):.0f}] m")

    # CPT
    dyn_cpt = str(OUTPUT_DIR / "_tmp_dyn_panels.cpt")
    dyn_range = (-float(vlim), float(vlim),
                 max(50.0, float(vlim) / 30.0))
    pygmt.makecpt(cmap=DYNTOPO_CPT, series=dyn_range, continuous=True,
                  background=True, output=dyn_cpt)

    fig = pygmt.Figure()
    pygmt.config(MAP_FRAME_TYPE="plain",
                 FONT_ANNOT_PRIMARY="8p,Helvetica,black",
                 FONT_LABEL="9p,Helvetica,black",
                 MAP_TITLE_OFFSET="-2p")

    # Winkel-Tripel aspect ratio ~ 1.637:1 (width:height).
    panel_h = width_cm / 1.637
    h_gap = 1.0        # vertical gap between rows
    v_gap = 0.4        # horizontal gap between cols

    # Origin starts at bottom-left.  Draw bottom row first; subsequent
    # shifts ADD to y as we go up.  The bottom row holds the younger
    # three ages (panels[3], panels[4], panels[5]).

    # Reserve room at the bottom for the colorbar.
    fig.shift_origin(xshift="0c", yshift="3.5c")

    # ---- BOTTOM ROW: ages[3], ages[4], ages[5] (e.g. 200, 100, 50 Ma) ----
    for col, idx in enumerate([3, 4, 5]):
        if col > 0:
            fig.shift_origin(xshift=f"{width_cm+v_gap}c")
        a, z, lat, lon = panels[idx]
        draw_panel(fig, to_xr(z, lat, lon), dyn_cpt,
                   proj=proj, age_label=f"{a} Ma")

    # Shared colorbar centered under the bottom row.
    cb_total_w = 3 * width_cm + 2 * v_gap
    fig.colorbar(projection=proj, region=REGION, cmap=dyn_cpt,
                 frame=['x+l"dynamic topography difference (m)"', "af"],
                 position=f"jBC+w{cb_total_w*0.55:.1f}c/0.3c"
                          f"+o-{width_cm + v_gap:.1f}c/-1.6c+h+e")

    # Back to the left column, shift up for the top row.
    fig.shift_origin(xshift=f"-{2*(width_cm+v_gap)}c",
                     yshift=f"{panel_h+h_gap}c")

    # ---- TOP ROW: ages[0], ages[1], ages[2] (e.g. 500, 400, 300 Ma) ----
    for col, idx in enumerate([0, 1, 2]):
        if col > 0:
            fig.shift_origin(xshift=f"{width_cm+v_gap}c")
        a, z, lat, lon = panels[idx]
        draw_panel(fig, to_xr(z, lat, lon), dyn_cpt,
                   proj=proj, age_label=f"{a} Ma")

    # ---- Save ----
    if write_png:
        out_png = OUTPUT_DIR / f"{out_basename}.png"
        if out_png.exists(): out_png.unlink()
        fig.savefig(out_png, dpi=dpi)
        print(f"  wrote {out_png}")
    if write_pdf:
        out_pdf = OUTPUT_DIR / f"{out_basename}.pdf"
        if out_pdf.exists(): out_pdf.unlink()
        fig.savefig(out_pdf)
        print(f"  wrote {out_pdf}")

    try: Path(dyn_cpt).unlink()
    except FileNotFoundError: pass


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--combined-dir", type=Path, default=DEFAULT_COMBINED_DIR,
                   help=f"directory of "
                        f"<age>Ma_corrected_plus_dyntopo_diff_young_SW.nc files "
                        f"(default: {DEFAULT_COMBINED_DIR})")
    p.add_argument("--ages", type=int, nargs=6, default=DEFAULT_AGES,
                   help=f"exactly 6 ages, top-row-first "
                        f"(default: {' '.join(str(a) for a in DEFAULT_AGES)})")
    p.add_argument("--vlim", type=float, default=DEFAULT_VLIM,
                   help="symmetric +/- range for the polar cpt; default 500 m")
    p.add_argument("--width-cm", type=float, default=7.0,
                   help="width of each panel in cm (default 7)")
    p.add_argument("--no-png", action="store_true")
    p.add_argument("--no-pdf", action="store_true")
    p.add_argument("--dpi", type=int, default=200)
    args = p.parse_args()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    out_basename = (f"dyntopo_diff_panels_"
                    f"{args.ages[0]}-{args.ages[-1]}Ma")
    make_figure(args.ages, args.combined_dir,
                out_basename, vlim=args.vlim,
                width_cm=args.width_cm,
                write_png=not args.no_png,
                write_pdf=not args.no_pdf,
                dpi=args.dpi)


if __name__ == "__main__":
    main()
