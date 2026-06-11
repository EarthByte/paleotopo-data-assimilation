#!/usr/bin/env python3
"""
=============================================================================
make_comparison_figures_dyntopo.py  —  Multi-panel pyGMT figures showing the
geochem-corrected S&W paleotopography next to the same maps with Young 2022
dynamic topography ADDITIVELY composed, plus the standalone dyntopo field.
=============================================================================

WHAT THIS DOES
    For each requested pair of ages, produces a 2-column x 3-row pyGMT
    figure in Winkel-Tripel projection.  Columns = ages, rows =
        row 1 - geochemically-corrected S&W elevation        (M_corrected)
        row 2 - corrected + Young 2022 dyntopo PER-STEP diff (M_combined)
        row 3 - Young 2022 dyntopo PER-STEP increment,       (z_dyntopo_diff
                dyntopo(t) - dyntopo(t - Δt), evaluated in         = correction signal
                plate frame then cookie-cut + rotated         in paleomag frame)
                to Scotese paleomag frame
    Shared colorbars at the bottom (one for the elevation rows, one for
    the dyntopo-difference row).

    The composition follows the methodology of
    `build_dyntopo_diff_correction.py`: present-day observed topography
    already contains today's dynamic-topography contribution, so the
    paleo correction is the time-difference dyntopo(t) - dyntopo(t - Δt).
    The subtraction is performed in plate reference frame (each cell
    rigidly attached to its continent) — subtracting at the same lat/lon
    in mantle/paleomag frame would compare different parts of different
    continents across time and is physically meaningless.  The
    plate-frame difference is then cookie-cut by Scotese 2023 continental
    polygons and rotated to the time-t paleomag frame via
    gplately.Raster.reconstruct, so it lines up with M_corrected.

    Default age pairs (one figure per pair) - constrained to the
    50-300 Ma window over which the Young 2022 dyntopo composition
    test is applied (see Supp methodology for the temporal cap
    rationale).  The 0 Ma slice is deliberately omitted: at t=0 the
    dyntopo correction is identically zero by construction, the
    "Corrected" label on row 1 does not apply (the present-day map is
    observed topography, not corrected), and the row 2 / row 3 panels
    convey no information.
        figure 1: (300 Ma, 250 Ma)
        figure 2: (200 Ma, 150 Ma)
        figure 3: (100 Ma,  50 Ma)

OUTPUT
    paths_scotese.OUTPUT_DIR / "SW_comparison_dyntopo_<age1>-<age2>Ma.png"
    paths_scotese.OUTPUT_DIR / "SW_comparison_dyntopo_<age1>-<age2>Ma.pdf"

USAGE
    cd <project>/scripts_Scotese
    python add_dyntopo_to_corrected_scotese.py --dyntopo-dir <...>   # produce inputs
    python make_comparison_figures_dyntopo.py                         # default 6 ages
    python make_comparison_figures_dyntopo.py --pairs 250 0           # 250 & 0
    python make_comparison_figures_dyntopo.py --no-pdf                # PNG only

OPTIONS
    --pairs  AGE1 AGE2 [AGE3 AGE4 ...]
        Explicit age list - even number of integers; consecutive pairs become
        figures.  Default: 250 200 150 100 50 0.
    --combined-dir PATH
        Directory of <age>Ma_corrected_plus_dyntopo_SW.nc files
        (default: PROJECT_ROOT/data/corrected_Scotese_plus_dyntopo).
    --no-pdf            skip PDF export
    --no-png            skip PNG export
    --dpi  INT          PNG dpi (default 200)
    --width-cm FLOAT    width of each map panel (default 10 cm)
    --vlim-dyntopo M    symmetric +/- range for the dyntopo cpt (default 1500)

DEPENDENCIES
    GMT 6.x, pygmt, xarray, netCDF4
    Optional: pygplates (only for plate-boundary overlay where resolvable -
                         0-100 Ma in S&W).  Set DRAW_SZ=False to skip.

NOTES
    - Inputs come from `add_dyntopo_to_corrected_scotese.py`.  That script
      must have been run first against your local Young 2022 / Scotese-frame
      dyntopo directory; the figure builder will refuse to run if the
      combined NetCDFs aren't present at the requested ages.
    - cpt range for the elevation rows is FIXED at +/-4000 m so the corrected
      and combined maps are directly comparable across panels (and across
      this figure vs the default Fig 5 comparison).
    - cpt range for the dyntopo row is FIXED at +/-`vlim-dyntopo` (default
      1500 m) using GMT's `polar` diverging cpt.
    - This script intentionally calls pyGMT for **all** plotting so the
      output is publication-grade Winkel-Tripel, matching the videos.
    - Scope-of-validity caveat for the caption: additive composition of
      geochem-corrected paleotopography with Young 2022 dyntopo will
      over-flood several continental interiors at multiple ages relative
      to the Kocsis & Scotese (ESR, 2021) paleocoastline compilation.
      This figure is intended as a diagnostic / supp-mat sanity check,
      NOT as the recommended composite paleotopography.
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
from paths_scotese import OUTPUT_DIR, PROJECT_ROOT


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
ELEV_CPT = "earth"
DYNTOPO_CPT = "polar"
ELEV_RANGE = (-4000.0, 4000.0, 250.0)      # series for makecpt (rows 1 & 2)
DEFAULT_VLIM_DYNTOPO = 1500.0              # +/- m for the dyntopo row (row 3)
PROJ_BASE = "R0"                            # Winkel-Tripel centred at lon=0
REGION = "g"                                # global
DRAW_SZ = True                              # overlay plate boundaries where resolvable

# Default age pairs constrained to 50-300 Ma (Young 2022 + S&W overlap
# minus the trivial t=0 slice).  Consecutive pairs of the list become
# figures: (300, 250), (200, 150), (100, 50).
DEFAULT_AGES = [300, 250, 200, 150, 100, 50]

COMBINED_FNAME_FMT = "{age}Ma_corrected_plus_dyntopo_diff_young_SW.nc"
DEFAULT_COMBINED_DIR = PROJECT_ROOT / "data" / "corrected_Scotese_plus_dyntopo_diff_young"


# ---------------------------------------------------------------------------
# Data helpers
# ---------------------------------------------------------------------------
def load_grid(t: int, combined_dir: Path):
    """Return (M_corrected, M_combined, z_dyntopo_diff, lat, lon, continent_mask).

    Reads the NetCDFs written by `build_dyntopo_diff_correction.py`:
    `M_corrected` (geochem-corrected S&W), `z_dyntopo_diff` (the
    plate-frame-computed time-difference dyntopo(t) - dyntopo(t - Δt)
    cookie-cut + rotated into Scotese paleomag frame at time t — i.e.
    the actual correction signal added to M_corrected), and
    `M_combined = M_corrected + z_dyntopo_diff` (with the land-guard
    clip applied where the correction pushed M_orig>0 cells below 0).

    `continent_mask` is NOT carried in the new NetCDFs — it's derived
    here from the original corrected-S&W NetCDF at the same age.
    """
    f = combined_dir / COMBINED_FNAME_FMT.format(age=int(t))
    if not f.exists():
        raise FileNotFoundError(
            f"{f} not found - run `python build_dyntopo_diff_correction.py "
            f"--source young --ages {int(t)}` first")
    with nc.Dataset(f) as d:
        Mc    = d.variables["M_corrected"][:].astype(float)
        Mcomb = d.variables["M_combined"][:].astype(float)
        z_dt  = d.variables["z_dyntopo_diff"][:].astype(float)
        lat   = d.variables["lat"][:].astype(float)
        lon   = d.variables["lon"][:].astype(float)
    # The build_dyntopo_diff_correction.py NetCDFs don't carry a
    # continent mask.  Derive it from the corresponding corrected-S&W
    # NetCDF (same age, no dyntopo involved) so plate-polygon coastline
    # contours still draw correctly.
    corrected_sw = (PROJECT_ROOT / "data" / "corrected_Scotese"
                    / f"{int(t)}Ma_corrected_SW.nc")
    if corrected_sw.exists():
        with nc.Dataset(corrected_sw) as d:
            cont = d.variables["continent_mask"][:].astype(bool)
    else:
        # Fallback: derive from M_corrected >= 0 (rough but usable)
        cont = (Mc >= 0)
    return Mc, Mcomb, z_dt, lat, lon, cont


def to_xr(arr: np.ndarray, lat, lon) -> xr.DataArray:
    """Wrap a (lat, lon) ndarray as an xarray DataArray, tagged as
    pixel-registered + geographic so GMT applies the correct periodic
    dateline wrap and doesn't emit the "Longitude range too small;
    geographic boundary condition changed to natural" warning.

    Mirrors the helper in make_comparison_figures.py.  Same justification:
    the S&W corrected NetCDFs are pixel-registered (cell-centre lons run
    -179.5..+179.5, centre-to-centre span = 359 deg even though the
    underlying coverage is 360 deg).
    """
    da = xr.DataArray(arr.astype(np.float32),
                      coords={"lat": lat, "lon": lon},
                      dims=("lat", "lon"), name="z")
    da.gmt.registration = 1   # 0 = gridline, 1 = pixel
    da.gmt.gtype = 1          # 0 = cartesian, 1 = geographic
    return da


def subduction_zones(t):
    """Return list of (lat, lon) arrays - empty if topologies don't
    resolve at t (S&W resolves topology only for ages <= ~100 Ma)."""
    if not DRAW_SZ:
        return []
    try:
        from plate_model_utils_scotese import subduction_zones as _sz
        return _sz(float(t))
    except Exception as e:
        print(f"  subduction zone overlay skipped at {t} Ma: {e}")
        return []


def plate_polygons(t):
    """Return list of closed (lat, lon) rings - one per resolved plate
    polygon - for the grey closed-boundary backstop on the figure
    panels.  Empty at ages where the S&W topology files don't resolve.
    Mirrors the same helper in make_comparison_figures.py.
    """
    if not DRAW_SZ:
        return []
    try:
        from plate_model_utils_scotese import plate_polygons as _pp
        return _pp(float(t))
    except Exception as e:
        print(f"  plate polygon overlay skipped at {t} Ma: {e}")
        return []


# ---------------------------------------------------------------------------
# Single panel drawer (one map within the multi-panel figure)
# ---------------------------------------------------------------------------
def draw_panel(fig, ds_da: xr.DataArray, cpt_path: str, title: str,
               proj: str, region: str = REGION,
               cont_da: xr.DataArray = None,
               sz_lines = None,
               plate_polys = None,
               show_age_label: str = None,
               show_coast: bool = True):
    """Render one map.  Caller is responsible for shifting origin.

    show_coast : draw the gray coastline contour derived from `cont_da`.
                 Switched off for the elevation rows where the terrain
                 colour ramp already makes the coastline visually obvious;
                 left on for the dyntopo row where the diverging colour
                 scheme on a mostly-near-zero field benefits from an
                 explicit continental outline.
    """
    # Note: title passed as `+t<title>` WITHOUT literal quotes - adding inner
    # double quotes (`+t"..."`) causes them to render visibly in the output.
    fig.basemap(region=region, projection=proj, frame=["af", f"+t{title}"])
    fig.grdimage(grid=ds_da, projection=proj, region=region,
                 cmap=cpt_path, nan_transparent=True)
    if show_coast and cont_da is not None and cont_da.values.any():
        fig.grdcontour(grid=cont_da, projection=proj, region=region,
                       levels=[0.5], pen="0.4p,gray30")
    # All plate boundaries as thin red lines on the figure maps - drawn
    # from the closed plate-polygon rings (resolved topologies), which
    # guarantees closure regardless of per-sub-segment feature-type
    # labelling.  The SZ overlay below is drawn ON TOP with a slightly
    # heavier pen so subduction zones still stand out from MORs /
    # transforms.  Matches the default make_comparison_figures.py style.
    if plate_polys:
        for ring in plate_polys:
            fig.plot(x=ring[:, 1], y=ring[:, 0],
                     pen="0.3p,red",
                     projection=proj, region=region)
    if sz_lines:
        for line in sz_lines:
            xy = np.column_stack([line[:, 1], line[:, 0]])
            fig.plot(x=xy[:, 0], y=xy[:, 1],
                     pen="0.6p,red", projection=proj, region=region)
    if show_age_label:
        fig.text(x=-170, y=65, text=show_age_label,
                 offset="-1.8c/0.8c", # 1.8 cm left, 0.8 cm up from anchor
                 font="14p,Helvetica-Bold,black",
                 region=region, projection=proj, justify="ML",
                 no_clip=True)


# ---------------------------------------------------------------------------
# Build one 2 x 3 comparison figure
# ---------------------------------------------------------------------------
def make_figure(age_left: int, age_right: int, out_basename: str,
                combined_dir: Path,
                vlim_dyntopo: float = DEFAULT_VLIM_DYNTOPO,
                width_cm: float = 10.0, write_png: bool = True,
                write_pdf: bool = True, dpi: int = 200):
    print(f"\n=== building dyntopo comparison figure: "
          f"{age_left} Ma | {age_right} Ma ===")
    proj = f"{PROJ_BASE}/{width_cm}c"  # e.g. "R0/10c"

    # Load both ages
    Mc_L, Mcomb_L, z_L, lat,  lon,  cont_L = load_grid(age_left,  combined_dir)
    Mc_R, Mcomb_R, z_R, lat2, lon2, cont_R = load_grid(age_right, combined_dir)
    assert lat.shape == lat2.shape and lon.shape == lon2.shape, "grid mismatch"

    sz_L = subduction_zones(age_left)
    sz_R = subduction_zones(age_right)
    pp_L = plate_polygons(age_left)
    pp_R = plate_polygons(age_right)

    # Coastline contour rasters derived from the rendered field for each
    # row.  Avoids any polygon vs DEM paleo-coastline mismatch.
    coast_L_corr_da = to_xr((Mc_L    >= 0).astype(np.float32), lat, lon)
    coast_R_corr_da = to_xr((Mc_R    >= 0).astype(np.float32), lat, lon)
    coast_L_comb_da = to_xr((Mcomb_L >= 0).astype(np.float32), lat, lon)
    coast_R_comb_da = to_xr((Mcomb_R >= 0).astype(np.float32), lat, lon)

    # Build the two CPT files
    elev_cpt = str(OUTPUT_DIR / "_tmp_elev_dyn.cpt")
    dyn_cpt  = str(OUTPUT_DIR / "_tmp_dyn.cpt")
    dyn_range = (-float(vlim_dyntopo), float(vlim_dyntopo),
                 max(50.0, float(vlim_dyntopo) / 30.0))
    pygmt.makecpt(cmap=ELEV_CPT,    series=ELEV_RANGE, continuous=True,
                  background=True, output=elev_cpt)
    pygmt.makecpt(cmap=DYNTOPO_CPT, series=dyn_range,  continuous=True,
                  background=True, output=dyn_cpt)

    fig = pygmt.Figure()
    pygmt.config(MAP_FRAME_TYPE="plain",
                 FONT_TITLE="12p,Helvetica-Bold,black",
                 FONT_ANNOT_PRIMARY="8p,Helvetica,black",
                 FONT_LABEL="9p,Helvetica,black",
                 MAP_TITLE_OFFSET="-2p")   # ~1.4 mm tighter than default 2p

    # Winkel-Tripel aspect ratio ~1.637:1 (width:height).
    panel_h = width_cm / 1.637
    h_gap = 1.4        # vertical gap between rows (in cm)
    v_gap = 0.6        # horizontal gap between cols (in cm)
    # Extra slot between rows 2 and 3 to host the shared elevation colorbar.
    cb_band = 1.4      # cm; height of the elevation-colorbar band

    # Origin starts at bottom-left; draw bottom row first so subsequent
    # shifts ADD to y as we go up.

    # ---- ROW 3 (bottom): dyntopo only ----
    # ASCII-only title (Ghostscript 10.x ISOLatin1+_Encoding workaround).
    fig.shift_origin(xshift="0c", yshift="3.5c")    # leave room for colorbars
    draw_panel(fig, to_xr(z_L, lat, lon), dyn_cpt,
               title="dyntopo diff (t - (t-Δt), plate frame)",
               proj=proj, cont_da=coast_L_corr_da, sz_lines=sz_L,
               plate_polys=pp_L,
               show_age_label=f"{age_left} Ma")
    fig.shift_origin(xshift=f"{width_cm+v_gap}c")
    draw_panel(fig, to_xr(z_R, lat, lon), dyn_cpt,
               title="dyntopo diff (t - (t-Δt), plate frame)",
               proj=proj, cont_da=coast_R_corr_da, sz_lines=sz_R,
               plate_polys=pp_R,
               show_age_label=f"{age_right} Ma")
    # dyntopo-anomaly colorbar across both bottom panels
    cb_total_w = 2*width_cm + v_gap
    fig.colorbar(projection=proj, region=REGION, cmap=dyn_cpt,
                 frame=['x+l"dynamic topography per-step increment (m)"', "af"],
                 position=f"jBC+w{cb_total_w*0.6:.1f}c/0.3c"
                          f"+o-{(width_cm+v_gap)/2:.1f}c/-1.6c+h+e")

    # Back to the left column for the middle row; y-shift bigger than the
    # normal inter-row gap so there's a `cb_band`-tall slot between rows
    # 2 and 3 hosting the shared elevation colorbar.
    fig.shift_origin(xshift=f"-{width_cm+v_gap}c",
                     yshift=f"{panel_h+h_gap+cb_band}c")

    # ---- ROW 2 (middle): corrected + dyntopo ----
    draw_panel(fig, to_xr(Mcomb_L, lat, lon), elev_cpt,
               title="Corrected + dyntopo diff",
               proj=proj, cont_da=coast_L_comb_da, sz_lines=sz_L,
               plate_polys=pp_L,
               show_age_label=f"{age_left} Ma",
               show_coast=False)
    fig.shift_origin(xshift=f"{width_cm+v_gap}c")
    draw_panel(fig, to_xr(Mcomb_R, lat, lon), elev_cpt,
               title="Corrected + dyntopo diff",
               proj=proj, cont_da=coast_R_comb_da, sz_lines=sz_R,
               plate_polys=pp_R,
               show_age_label=f"{age_right} Ma",
               show_coast=False)
    # elevation colorbar in the cb_band gap between rows 2 and 3
    fig.colorbar(projection=proj, region=REGION, cmap=elev_cpt,
                 frame=['x+l"elevation (m)"', "af"],
                 position=f"jBC+w{cb_total_w*0.6:.1f}c/0.3c"
                          f"+o-{(width_cm+v_gap)/2:.1f}c/-1.4c+h+e")

    # ---- ROW 1 (top): corrected ----
    fig.shift_origin(xshift=f"-{width_cm+v_gap}c",
                     yshift=f"{panel_h+h_gap}c")
    draw_panel(fig, to_xr(Mc_L, lat, lon), elev_cpt,
               title="Corrected",
               proj=proj, cont_da=coast_L_corr_da, sz_lines=sz_L,
               plate_polys=pp_L,
               show_age_label=f"{age_left} Ma",
               show_coast=False)
    fig.shift_origin(xshift=f"{width_cm+v_gap}c")
    draw_panel(fig, to_xr(Mc_R, lat, lon), elev_cpt,
               title="Corrected",
               proj=proj, cont_da=coast_R_corr_da, sz_lines=sz_R,
               plate_polys=pp_R,
               show_age_label=f"{age_right} Ma",
               show_coast=False)

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

    # Clean up temporary CPTs
    for f in [elev_cpt, dyn_cpt]:
        try: Path(f).unlink()
        except FileNotFoundError: pass


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pairs", type=int, nargs="+", default=DEFAULT_AGES,
                   help="ages to compare, in pairs (e.g. 250 200 150 100 50 0)")
    p.add_argument("--combined-dir", type=Path, default=DEFAULT_COMBINED_DIR,
                   help="directory of <age>Ma_corrected_plus_dyntopo_SW.nc files")
    p.add_argument("--no-png", action="store_true")
    p.add_argument("--no-pdf", action="store_true")
    p.add_argument("--dpi", type=int, default=200)
    p.add_argument("--width-cm", type=float, default=10.0,
                   help="width of each map panel in cm (default 10)")
    p.add_argument("--vlim-dyntopo", type=float, default=DEFAULT_VLIM_DYNTOPO,
                   help="symmetric +/- range for the dyntopo cpt; default 1500 m")
    args = p.parse_args()

    ages = args.pairs
    if len(ages) % 2 != 0:
        raise SystemExit("ERROR: --pairs must contain an even number of ages")

    # Validate that requested ages have combined NetCDFs
    combined_dir = args.combined_dir
    missing = [a for a in ages
               if not (combined_dir / COMBINED_FNAME_FMT.format(age=a)).exists()]
    if missing:
        raise SystemExit(
            f"ERROR: no combined NetCDF for ages {missing} in {combined_dir}. "
            f"Run `python add_dyntopo_to_corrected_scotese.py --dyntopo-dir <...>` "
            f"first.")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for i in range(0, len(ages), 2):
        a_l, a_r = ages[i], ages[i + 1]
        out_basename = f"SW_comparison_dyntopo_{a_l}-{a_r}Ma"
        make_figure(a_l, a_r, out_basename, combined_dir,
                    vlim_dyntopo=args.vlim_dyntopo,
                    width_cm=args.width_cm,
                    write_png=not args.no_png,
                    write_pdf=not args.no_pdf,
                    dpi=args.dpi)


if __name__ == "__main__":
    main()
