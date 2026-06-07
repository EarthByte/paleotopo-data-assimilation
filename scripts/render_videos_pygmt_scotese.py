#!/usr/bin/env python3
"""
=============================================================================
render_videos_pygmt_scotese.py — Publication videos (Winkel-Tripel) for S&W
=============================================================================

Renders up to four MP4 videos from the Scotese & Wright corrected NetCDFs
using pyGMT in Winkel-Tripel projection (GMT 'globe' and 'polar' CPTs).
Publication-grade counterpart of render_videos_cartopy_scotese.py.

OUTPUT VIDEOS  (paths_scotese.OUTPUT_DIR)
    SW_paleotopo_elevation_<...>Ma.mp4           corrected elevation
    SW_paleotopo_original_<...>Ma.mp4            original S&W input
    SW_paleotopo_delta_<...>Ma.mp4               Δz = corrected − input
    SW_paleotopo_elevation_samples_<...>Ma.mp4   corrected + geochem controls

PER-FRAME LAYOUT
    - Winkel-Tripel projection (R0/18c)
    - elevation cpt: globe, range −4000..+4000 m (matches Fig 5 comparison)
    - delta cpt:     polar, range −2000..+2000 m
    - plate-boundary topology overlay in DELTA mode only, drawn in BLACK
      (matches the Merdith delta video).  The S&W plate model only
      resolves topologies for ages ≲ 100 Ma — `plate_model_utils_scotese
      .topology_lines(t)` returns empty for older slices and the
      overlay silently degrades to nothing.  Elevation / original /
      elevation_samples modes carry no overlay.
    - age annotation TOP-LEFT (matches Merdith convention)
    - 200 dpi PNG → libx264 yuv420p MP4

INPUT
    paths_scotese.CORRECTED_DIR / "<age>Ma_corrected_SW.nc"

DEPENDENCIES
    GMT 6.x, pygmt, xarray, netCDF4, ffmpeg, pygplates (optional)

USAGE
    cd <project>/scripts
    python render_videos_pygmt_scotese.py                      # all four, all ages
    python render_videos_pygmt_scotese.py --fps 12
    python render_videos_pygmt_scotese.py --modes original elevation_samples
    python render_videos_pygmt_scotese.py --ages 0 50 100 200 500

OPTIONS  (same flags as the cartopy version)
    --cadence INT     render every N Ma (default: every S&W age)
    --ages   INT ...  explicit age list
    --fps    INT      video frame rate (default 10)
    --modes  ...      which video(s) to build
    --keep-frames     don't delete the per-frame PNG directory after stitching

RUNTIME (1° / 200 dpi)
    ~3 s per frame.  Full 107-frame × 4-mode sweep ≈ 20 min per video.
=============================================================================
"""
from __future__ import annotations
import os, sys, argparse, subprocess, shutil
from pathlib import Path
import numpy as np
import netCDF4 as nc

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from paths_scotese import CORRECTED_DIR, OUTPUT_DIR, CSV_PATH
from sw_io import available_ages as sw_available_ages

PROJ_ROOT = HERE.parent
FRAME_DIR = OUTPUT_DIR / "video_frames_pygmt_SW"
OUT_DIR = OUTPUT_DIR
FRAME_DIR.mkdir(parents=True, exist_ok=True)
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Plate-boundary topology overlay is drawn ONLY on the delta video,
# in black, to match the Merdith delta convention.  Setting False
# disables the overlay completely (e.g. if pygplates isn't installed).
DRAW_SZ = True

import pygmt
import xarray as xr


# ---------------------------------------------------------------------------
# Adaptive temporal half-window — must match assimilate_scotese.dt_half()
# ---------------------------------------------------------------------------
def _dt_half(t):
    if t < 200:  return 5.0
    if t < 500:  return 10.0
    if t < 800:  return 20.0
    return 30.0


# Geochem cache + per-frame reconstruction (S&W plate IDs via continental
# polygons; oceanic samples dropped — see sample_reconstruct_scotese.py).
_GEOCHEM_CACHE = None
_SAMPLE_RECONSTRUCTOR = None
_PLATE_ID_CACHE = CSV_PATH.parent / "sample_plate_ids_SW.npy"


def _geochem():
    global _GEOCHEM_CACHE, _SAMPLE_RECONSTRUCTOR
    if _GEOCHEM_CACHE is None:
        import pandas as pd
        from sample_reconstruct_scotese import ScoteseSampleReconstructor
        if _SAMPLE_RECONSTRUCTOR is None:
            _SAMPLE_RECONSTRUCTOR = ScoteseSampleReconstructor()
        df = pd.read_csv(CSV_PATH, low_memory=False)
        mid_cols = [
            "Isostatic_Elevation_absolute_km",
            "Brown_Isostatic_Elevation_absolute_km",
            "Davis_Isostatic_Elevation_absolute_km",
            "Condie_Isostatic_Elevation_absolute_km",
            "Herz_0.08_Isostatic_Elevation_absolute_km",
            "Herz_0.38 Isostatic_Elevation_absolute_km",
        ]
        df = df.dropna(subset=["Lat", "Lon", "Age_Ma"])
        df = df[df.get("Missing_Risk_DeltaT (km)", 0).fillna(0) <= 5].copy()
        df["z_obs_m"] = df[mid_cols].astype(float).mean(axis=1, skipna=True) * 1000.0
        df = df.dropna(subset=["z_obs_m"]).reset_index(drop=True)
        if _PLATE_ID_CACHE.exists():
            pids = np.load(_PLATE_ID_CACHE)
            if len(pids) == len(df):
                df["sw_plate_id"] = pids
        if "sw_plate_id" not in df.columns:
            df["sw_plate_id"] = _SAMPLE_RECONSTRUCTOR.assign_plate_ids(df)
            np.save(_PLATE_ID_CACHE, df["sw_plate_id"].to_numpy())
        df = df[df["sw_plate_id"] != 0].reset_index(drop=True)
        _GEOCHEM_CACHE = df[["Lat", "Lon", "Age_Ma", "z_obs_m", "sw_plate_id"]]
    return _GEOCHEM_CACHE


def _samples_at_age(t, half, lat=None, lon=None, cont=None,
                    M=None, Mc=None):
    """Samples in the ±half window reconstructed to age t.

    If lat/lon/cont are provided, drop samples whose reconstructed
    position lands outside the per-slice continent mask (terranes whose
    paleo-rotation has broken down at age t — polygon valid-time
    mismatch, missing terrane rotation, etc.).

    If M (original) and/or Mc (corrected) elevation arrays are also
    provided, drop samples landing on flooded cells (elevation ≤ 0 m).
    This catches the rare case of a sample reconstructing into a
    drowned-interior basin or flooded shelf cell that is inside the
    polygon footprint but would visually plot over the ocean in the
    rendered map."""
    from sw_io import nearest_cell_index
    df = _geochem()
    sub = df[(df["Age_Ma"] >= t - half) & (df["Age_Ma"] <= t + half)].reset_index(drop=True)
    if sub.empty:
        return sub.assign(rlat_t=[], rlon_t=[])
    rlat_t, rlon_t = _SAMPLE_RECONSTRUCTOR.reconstruct(sub, t, plate_col="sw_plate_id")
    sub = sub.copy()
    sub["rlat_t"] = rlat_t
    sub["rlon_t"] = rlon_t
    sub = sub.dropna(subset=["rlat_t", "rlon_t"]).reset_index(drop=True)
    if cont is not None and lat is not None and lon is not None and len(sub):
        iy = nearest_cell_index(lat, sub["rlat_t"].values)
        ix = nearest_cell_index(lon, sub["rlon_t"].values)
        keep = cont[iy, ix].copy()
        if M is not None:
            keep &= (M[iy, ix] > 0)
        if Mc is not None:
            keep &= (Mc[iy, ix] > 0)
        sub = sub[keep].reset_index(drop=True)
    return sub


def load_grid_xr(t: int) -> xr.Dataset:
    f = CORRECTED_DIR / f"{int(t)}Ma_corrected_SW.nc"
    return xr.open_dataset(f)


def make_xr_dataarray(arr: np.ndarray, lat, lon, name="z") -> xr.DataArray:
    return xr.DataArray(arr, coords={"lat": lat, "lon": lon},
                        dims=("lat", "lon"), name=name)


# Elevation modes use a FIXED range of -4000..+4000 m so the videos can
# be compared frame-to-frame and slice-to-slice with the same colour
# mapping.  Data outside the range overflows into the colorbar
# end-triangles enabled via `+e` in the colorbar `position` argument
# below.  Matches the comparison-figure cpt range.
ELEVATION_RANGE = (-4000.0, 4000.0)
DATA_RANGE = {
    "elevation": ELEVATION_RANGE,
    "elevation_samples": ELEVATION_RANGE,
    "original": ELEVATION_RANGE,
    "delta": (-2000.0, 2000.0),
}


def _nice_step(span: float) -> float:
    """Pick a colour-table step that yields ~25-50 intervals across the span."""
    for s in (50, 100, 200, 250, 500, 1000):
        if span / s <= 50:
            return s
    return 2000.0


def compute_data_range(ages, modes):
    """Only computes the delta range from data — elevation is fixed."""
    if "delta" not in modes:
        return
    mn_d, mx_d = np.inf, -np.inf
    for t in ages:
        with nc.Dataset(CORRECTED_DIR / f"{int(t)}Ma_corrected_SW.nc") as d:
            a = d.variables["delta"][:]
            mn_d = min(mn_d, float(np.nanmin(a))); mx_d = max(mx_d, float(np.nanmax(a)))
    d_abs = max(abs(mn_d), abs(mx_d), 1.0)
    step = _nice_step(2 * d_abs)
    d_abs = float(np.ceil(d_abs / step) * step)
    DATA_RANGE["delta"] = (-d_abs, d_abs)
    print(f"  delta cpt range derived from data: ±{d_abs:.0f} m  "
          f"(elevation range fixed at {ELEVATION_RANGE[0]:.0f} … {ELEVATION_RANGE[1]:.0f} m)")


def render_frame(t: int, mode: str, frame_path: Path):
    """Render one frame. mode ∈ {'elevation', 'original', 'delta', 'elevation_samples'}."""
    ds = load_grid_xr(t)
    lat = ds["lat"].values
    lon = ds["lon"].values
    overlay_samples = False
    elev_lo, elev_hi = DATA_RANGE["elevation"]
    elev_step = _nice_step(elev_hi - elev_lo)
    d_lo, d_hi = DATA_RANGE["delta"]
    d_step = _nice_step(d_hi - d_lo)
    # Note: per-frame top titles were removed in 2026-05-26 — the only
    # remaining text annotation is the time stamp in the upper-left
    # corner of the map.  `cb_label` is still ASCII-only because
    # Ghostscript 10.x doesn't define the ISOLatin1+_Encoding GMT emits
    # for non-ASCII glyphs (Greek delta, true minus, plus-minus).
    if mode == "elevation":
        arr = ds["M_corrected"].values
        cpt = "earth"; series = (elev_lo, elev_hi, elev_step)
        cb_label = "elevation (m)"
    elif mode == "elevation_samples":
        arr = ds["M_corrected"].values
        cpt = "earth"; series = (elev_lo, elev_hi, elev_step)
        cb_label = "elevation (m)"
        overlay_samples = True
    elif mode == "original":
        arr = ds["M_orig"].values
        cpt = "earth"; series = (elev_lo, elev_hi, elev_step)
        cb_label = "elevation (m)"
    elif mode == "delta":
        # Mask Δz to the actual elevation-field paleo-land footprint
        # (M_orig ≥ 0 ∪ M_corrected ≥ 0), NOT the Scotese 2008
        # continent-polygon raster.  Otherwise cells that are inside
        # the polygon footprint but below sea level in the rendered
        # S&W18 topography paint a faint Δz "ghost continent" outside
        # the real paleo-coastline — particularly visible before
        # ~350 Ma where the polygon vs DEM mismatch is large.
        _Mo = ds["M_orig"].values
        _Mc = ds["M_corrected"].values
        _land = (_Mo >= 0) | (_Mc >= 0)
        arr = np.where(_land, ds["delta"].values, np.nan)
        cpt = "polar"; series = (d_lo, d_hi, d_step)
        cb_label = "delta z (m)"
        overlay_samples = True
    else:
        raise ValueError(mode)

    da = make_xr_dataarray(arr.astype(np.float32), lat, lon)
    # Tell GMT the grid is pixel-registered (cell-centre lons run
    # -179.5..+179.5, so centre-to-centre span is 359° even though the
    # underlying coverage is 360°) and geographic (so grdimage applies
    # periodic dateline wrap).  Without these two attrs GMT defaults
    # to gridline-registered Cartesian, sees the 359° span literally,
    # and emits the "Longitude range too small; geographic boundary
    # condition changed to natural" warning every frame.
    da.gmt.registration = 1   # 0 = gridline, 1 = pixel
    da.gmt.gtype = 1          # 0 = cartesian, 1 = geographic
    fig = pygmt.Figure()
    pygmt.config(MAP_FRAME_TYPE="plain",
                 FONT_ANNOT_PRIMARY="9p,Helvetica,black",
                 FONT_LABEL="10p,Helvetica,black",
                 FONT_TITLE="14p,Helvetica-Bold,black")
    proj = "R0/18c"
    region = "g"
    # No title — frame is just the axis annotation/tick scaffolding.
    fig.basemap(region=region, projection=proj, frame="af")
    pygmt.makecpt(cmap=cpt, series=series, continuous=True,
                  background=True, output="paleo.cpt")
    fig.grdimage(grid=da, projection=proj, region=region,
                 cmap="paleo.cpt", nan_transparent=True)

    # No coastline contour on the videos — at 1° / 5 Myr the contour
    # outlines pixel-edge artefacts that draw the eye away from the
    # actual elevation signal.  Static paper figures keep the contour
    # because they're read at a different scale.
    cont_mask = ds["continent_mask"].values.astype(bool)

    # Plate-boundary topology overlay in DELTA mode only, drawn in
    # BLACK to match the Merdith delta video.  The S&W plate model
    # only resolves topologies for ages ≲ 100 Ma — `topology_lines`
    # returns empty for older slices and the overlay silently
    # degrades to nothing.  Elevation modes carry no overlay.
    if DRAW_SZ and mode == "delta":
        try:
            from plate_model_utils_scotese import topology_lines
            from topology_render import draw_topologies_pygmt
        except ImportError as _e:
            if not getattr(render_frame, "_topology_warned", False):
                print(f"  topology overlay disabled: {_e}")
                render_frame._topology_warned = True
        else:
            draw_topologies_pygmt(
                fig, topology_lines(float(t)),
                projection=proj, region=region,
                other_pen="0.4p,black",
                sz_pen="0.7p,black",
                sz_fill="black",
                polygon_pen="0.25p,black",
            )

    if overlay_samples:
        # Samples are RECONSTRUCTED to age t (not at static depositional position).
        M_orig = ds["M_orig"].values
        M_corr = ds["M_corrected"].values
        sub = _samples_at_age(t, _dt_half(t), lat=lat, lon=lon,
                              cont=cont_mask, M=M_orig, Mc=M_corr)
        if len(sub):
            if mode == "delta":
                # Δz mode: position-only marker (open black circles).
                # Colour-coding by z_obs on a diverging Δz cmap would
                # mislead — z_obs is absolute, Δz is a signed correction.
                fig.plot(x=sub["rlon_t"].values, y=sub["rlat_t"].values,
                         style="c0.10c", pen="0.5p,black",
                         projection=proj, region=region)
            else:
                # Elevation mode: filled circles coloured by z_obs on the
                # same terrain cmap as the background.
                fig.plot(x=sub["rlon_t"].values, y=sub["rlat_t"].values,
                         fill=sub["z_obs_m"].values, cmap="paleo.cpt",
                         style="c0.10c", pen="0.3p,black",
                         projection=proj, region=region)
        fig.text(x=-170, y=80, text=f"n = {len(sub)} samples",
                 font="11p,Helvetica,black", region=region, projection=proj,
                 justify="ML", fill="white", pen="0.25p,gray50")

    # `+e` on the colorbar position draws triangle arrows at BOTH ends
    # indicating the data extends past the fixed ±4000 m cpt range.
    fig.colorbar(projection=proj, region=region, cmap="paleo.cpt",
                 frame=[f'x+l"{cb_label}"'], position="JBC+w12c/0.3c+o0/1c+h+e")
    # Age stamp in the upper-left corner — matches the Merdith pyGMT
    # convention.  `position="TL"` anchors at the top-left of the map
    # bounding box, `offset` nudges 0.15 cm LEFT and 0.5 cm DOWN (GMT
    # "+x is right, +y is up") so the stamp sits just inside the
    # upper-left of the map, below the frame border.  History:
    # +0.9c → +0.3c → -0.7c → -0.4c → -0.5c as the title was removed
    # and the stamp was tuned to its final position.  `no_clip=True`
    # lets the text draw outside the projection clip.
    fig.text(text=f"{int(t)} Ma",
             font="22p,Helvetica-Bold,black",
             position="TL", justify="TL",
             offset="-0.15c/-0.5c",
             no_clip=True,
             region=region, projection=proj)
    fig.savefig(frame_path, dpi=200)
    if Path("paleo.cpt").exists():
        Path("paleo.cpt").unlink()


def stitch_video(frame_pattern: Path, out_path: Path, fps: int, n_frames: int):
    """Explicit -frames:v n_frames so stale frames left by a previous
    cadence/age-range run don't get appended."""
    cmd = ["ffmpeg", "-y", "-framerate", str(fps),
           "-i", str(frame_pattern),
           "-frames:v", str(n_frames),
           "-c:v", "libx264", "-pix_fmt", "yuv420p",
           "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2",
           "-crf", "20", str(out_path)]
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--cadence", type=int, default=None)
    p.add_argument("--ages", type=int, nargs="*")
    p.add_argument("--fps", type=int, default=10)
    p.add_argument("--modes", nargs="+",
                   default=["elevation", "delta", "original", "elevation_samples"],
                   choices=["elevation", "delta", "original", "elevation_samples"])
    p.add_argument("--keep-frames", action="store_true")
    p.add_argument("--force", action="store_true",
                   help="wipe cached per-frame PNGs (under video_frames_pygmt_SW/<mode>/) "
                        "before rendering. Only matters when --keep-frames is set or a "
                        "previous run with --keep-frames left behind cached frames; "
                        "without --keep-frames, frames are deleted after MP4 stitching.")
    args = p.parse_args()

    if args.ages:
        ages = sorted(args.ages, reverse=True)
    else:
        all_ages = sw_available_ages()
        if args.cadence:
            ages = sorted([a for a in all_ages if a % args.cadence == 0], reverse=True)
        else:
            ages = sorted(all_ages, reverse=True)

    print(f"Scanning {len(ages)} NetCDFs to derive cpt range from actual data …")
    compute_data_range(ages, args.modes)

    for mode in args.modes:
        sub = FRAME_DIR / mode
        sub.mkdir(exist_ok=True)
        existing = list(sub.glob("frame_*.png"))
        if args.force and existing:
            print(f"[{mode}] --force: wiping {len(existing)} cached frame(s)")
            for f in existing:
                f.unlink()
        elif existing and len(existing) >= len(ages):
            print(f"[{mode}] reusing {len(existing)} cached frames "
                  f"(pass --force to re-render after a code change)")
        for n, t in enumerate(ages):
            frame = sub / f"frame_{n:04d}.png"
            if frame.exists(): continue
            print(f"[{mode}] rendering frame {n}/{len(ages)} at {t} Ma")
            render_frame(t, mode, frame)
        out_video = OUT_DIR / f"SW_paleotopo_{mode}_{ages[0]}-{ages[-1]}Ma.mp4"
        if out_video.exists(): out_video.unlink()
        stitch_video(sub / "frame_%04d.png", out_video, args.fps, n_frames=len(ages))
        print(f"  wrote {out_video}")
        if not args.keep_frames:
            shutil.rmtree(sub)


if __name__ == "__main__":
    main()
