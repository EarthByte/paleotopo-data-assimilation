#!/usr/bin/env python3
"""
=============================================================================
render_videos_pygmt_dyntopo_scotese.py  —  Publication videos (Winkel-Tripel)
                                            of the corrected + dyntopo
                                            paleotopography composition
=============================================================================

Renders publication-grade Winkel-Tripel MP4 videos from the combined NetCDFs
produced by build_dyntopo_diff_correction.py:

    M_combined(t)        = M_corrected(t) + Δz_paleomag(t)
    Δz_paleomag(t)       = z_dyntopo(t) − z_dyntopo(t − Δt)  rotated to paleomag

Two modes are produced by default:

    combined            M_combined — the geochem-corrected paleotopography
                        with the per-step dyntopo increment added (the
                        headline dyntopo-composition video).
    dyntopo_diff        z_dyntopo_diff alone — the per-step dyntopo
                        increment applied at each age, in Scotese paleomag
                        frame.

Two additional modes are available on request:

    corrected           M_corrected from the same combined NetCDFs, rendered
                        in the identical projection / colour scale /
                        age-stamp layout so paired videos can be compared
                        frame-by-frame against the combined output.
    dyntopo_absolute    z_dyntopo_paleomag — the absolute past dynamic
                        topography on the continents at each age (Young
                        2022 plate-frame field cookie-cut by Scotese 2023
                        polygons and rotated into the time-t paleomag
                        frame).  Polar cpt at ±1200 m.  NOT to be added
                        directly to M_corrected (that would double-count
                        today's dyntopo contribution); provided as the
                        absolute-state visualisation that accompanies the
                        per-step diff in M_combined.

PER-FRAME LAYOUT
    - Winkel-Tripel projection (R0/18c)
    - combined / corrected cpt: GMT 'earth', range −4000..+4000 m
      (matches the main-text Fig 5 and the default Scotese paleotopo videos)
    - dyntopo_diff cpt:         GMT 'polar', range ±500 m
    - dyntopo_absolute cpt:     GMT 'polar', range ±1200 m
    - plate-boundary overlay in BLACK on the dyntopo_diff and
      dyntopo_absolute videos
      (matches the delta-mode convention in render_videos_pygmt_scotese.py;
      S&W topologies resolve only for ages ≲ 100 Ma, older slices silently
      degrade to no overlay)
    - age annotation in the upper-left corner (matches Merdith / Scotese
      pyGMT convention)
    - no per-frame title (just the age stamp)
    - 200 dpi PNG → libx264 yuv420p MP4

INPUT
    PROJECT_ROOT / "data" / "corrected_Scotese_plus_dyntopo_diff_young"
        / "<age>Ma_corrected_plus_dyntopo_diff_young_SW.nc"
    Variables read: M_corrected, M_combined, z_dyntopo_diff,
                    z_dyntopo_paleomag.

OUTPUT VIDEOS  (paths_scotese.OUTPUT_DIR)
    SW_paleotopo_combined_<a1>-<a2>Ma.mp4
    SW_paleotopo_dyntopo_diff_<a1>-<a2>Ma.mp4
    SW_paleotopo_corrected_<a1>-<a2>Ma.mp4         (optional)
    SW_paleotopo_dyntopo_absolute_<a1>-<a2>Ma.mp4  (optional)

USAGE
    cd <project>/scripts_Scotese
    # Composition NetCDFs must exist first
    python build_dyntopo_diff_correction.py --source young \\
        --ages $(seq 0 5 300) --step-myr 5
    # Then render the videos
    python render_videos_pygmt_dyntopo_scotese.py
    python render_videos_pygmt_dyntopo_scotese.py --modes combined dyntopo_diff
    python render_videos_pygmt_dyntopo_scotese.py --fps 12 --ages 0 50 100 200 300

OPTIONS
    --combined-dir PATH    directory of
                           <age>Ma_corrected_plus_dyntopo_diff_young_SW.nc
                           (default: PROJECT_ROOT/data/corrected_Scotese_plus_dyntopo_diff_young)
    --ages INT ...         explicit list of ages (default: all available)
    --fps INT              video frame rate (default 10)
    --modes ...            which video(s) to build: combined / dyntopo_diff /
                           corrected (default: combined dyntopo_diff)
    --keep-frames          don't delete the per-frame PNG directory after stitching
    --force                wipe the per-mode frame cache before rendering

DEPENDENCIES
    GMT 6.x, pygmt, xarray, netCDF4, ffmpeg, pygplates (optional for overlay)
=============================================================================
"""
from __future__ import annotations
import argparse, os, re, shutil, subprocess, sys
from pathlib import Path
import numpy as np
import netCDF4 as nc

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from paths_scotese import OUTPUT_DIR, PROJECT_ROOT, VIDEOS_DIR

PROJ_ROOT = HERE.parent
DEFAULT_COMBINED_DIR = PROJECT_ROOT / "data" / "corrected_Scotese_plus_dyntopo_diff_young"
FRAME_DIR = OUTPUT_DIR / "video_frames_pygmt_SW_dyntopo"
OUT_DIR = OUTPUT_DIR
FRAME_DIR.mkdir(parents=True, exist_ok=True)
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Plate-boundary topology overlay is drawn ONLY on the dyntopo_diff video,
# in black (matching the Merdith / Scotese delta-mode convention).
DRAW_SZ = True

import pygmt
import xarray as xr

COMBINED_FNAME_FMT = "{age}Ma_corrected_plus_dyntopo_diff_young_SW.nc"
COMBINED_FNAME_RE  = re.compile(r"^(\d+)Ma_corrected_plus_dyntopo_diff_young_SW\.nc$")

# Globally-consistent cpt series.  Elevation modes use a FIXED range of
# −4000 to +4000 m to match render_videos_pygmt_scotese.py and the main-text
# Fig 5 comparison.  The dyntopo_diff (per-step) range is fixed at ±500 m
# and the dyntopo_absolute range at ±1200 m — the per-step magnitudes are
# tens of metres, while the absolute past dyntopo can reach a kilometre
# or so at major upwelling / subsidence regions, so ±1200 m gives good
# contrast without saturating the extremes.
ELEVATION_RANGE = (-4000.0, 4000.0)
DYNTOPO_DIFF_RANGE = (-500.0, 500.0)
DYNTOPO_ABS_RANGE  = (-1200.0, 1200.0)


def auto_discover_ages(combined_dir: Path) -> list[int]:
    if not combined_dir.is_dir():
        return []
    return sorted(
        int(COMBINED_FNAME_RE.match(p.name).group(1))
        for p in combined_dir.iterdir()
        if COMBINED_FNAME_RE.match(p.name)
    )


def load_grid_xr(combined_dir: Path, t: int, var: str) -> xr.DataArray:
    f = combined_dir / COMBINED_FNAME_FMT.format(age=int(t))
    if not f.exists():
        raise FileNotFoundError(
            f"{f} not found - run `python build_dyntopo_diff_correction.py "
            f"--source young --ages {int(t)} --step-myr 5` first")
    ds = xr.open_dataset(f)
    da = ds[var]
    return da


def _nice_step(span: float) -> float:
    """Pick a cpt step that yields ~25-50 intervals across the span."""
    for s in (50, 100, 200, 250, 500, 1000):
        if span / s <= 50:
            return s
    return 2000.0


def render_frame(combined_dir: Path, t: int, mode: str, frame_path: Path):
    """Render one frame.

    mode ∈ {'combined', 'corrected', 'dyntopo_diff', 'dyntopo_absolute'}.
    """
    if mode == "combined":
        var = "M_combined"
        cpt = "earth"
        lo, hi = ELEVATION_RANGE
        cb_label = "elevation (m)"
    elif mode == "corrected":
        var = "M_corrected"
        cpt = "earth"
        lo, hi = ELEVATION_RANGE
        cb_label = "elevation (m)"
    elif mode == "dyntopo_diff":
        var = "z_dyntopo_diff"
        cpt = "polar"
        lo, hi = DYNTOPO_DIFF_RANGE
        cb_label = "dyntopo per-step increment (m)"
    elif mode == "dyntopo_absolute":
        var = "z_dyntopo_paleomag"
        cpt = "polar"
        lo, hi = DYNTOPO_ABS_RANGE
        cb_label = "dynamic topography (m)"
    else:
        raise ValueError(mode)
    step = _nice_step(hi - lo)
    series = (lo, hi, step)

    da = load_grid_xr(combined_dir, t, var).astype(np.float32)
    # For both dyntopo modes, mask to the land footprint of M_combined /
    # M_corrected so the field doesn't paint over the open ocean where its
    # value is meaningless (cookie-cut is by continental polygons, but
    # cells just off the coast still carry small non-zero values from the
    # rotation step).
    if mode in ("dyntopo_diff", "dyntopo_absolute"):
        Mc = load_grid_xr(combined_dir, t, "M_corrected").values
        Mcomb = load_grid_xr(combined_dir, t, "M_combined").values
        land = (Mc >= 0) | (Mcomb >= 0)
        da = xr.DataArray(np.where(land, da.values, np.nan),
                          coords=da.coords, dims=da.dims, name=da.name)
    # Tag as pixel-registered + geographic so GMT applies periodic dateline
    # wrap and doesn't print the "Longitude range too small" warning.
    da.gmt.registration = 1
    da.gmt.gtype = 1

    fig = pygmt.Figure()
    pygmt.config(MAP_FRAME_TYPE="plain",
                 FONT_ANNOT_PRIMARY="9p,Helvetica,black",
                 FONT_LABEL="10p,Helvetica,black",
                 FONT_TITLE="14p,Helvetica-Bold,black")
    proj = "R0/18c"
    region = "g"
    fig.basemap(region=region, projection=proj, frame="af")
    pygmt.makecpt(cmap=cpt, series=series, continuous=True,
                  background=True, output="paleo.cpt")
    fig.grdimage(grid=da, projection=proj, region=region,
                 cmap="paleo.cpt", nan_transparent=True)

    # Plate-boundary overlay on dyntopo_diff only, in black (matches the
    # default Scotese delta-mode convention).  S&W topologies resolve only
    # for ages ≲ 100 Ma; older slices silently get no overlay.
    if DRAW_SZ and mode in ("dyntopo_diff", "dyntopo_absolute"):
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

    # `+e` colorbar end-triangles indicate cells outside the fixed range.
    fig.colorbar(projection=proj, region=region, cmap="paleo.cpt",
                 frame=[f'x+l"{cb_label}"'],
                 position="JBC+w12c/0.3c+o0/1c+h+e")
    # Age stamp in the upper-left corner of the map.  Matches the convention
    # in render_videos_pygmt_scotese.py: position="TL", offset 0.15 cm left
    # and 0.5 cm down from the top-left bbox corner.
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
    """Pin -frames:v so leftover frames from a previous run can't sneak in."""
    cmd = ["ffmpeg", "-y", "-framerate", str(fps),
           "-i", str(frame_pattern),
           "-frames:v", str(n_frames),
           "-c:v", "libx264", "-pix_fmt", "yuv420p",
           "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2",
           "-crf", "20", str(out_path)]
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--combined-dir", type=Path, default=DEFAULT_COMBINED_DIR,
                   help="directory of "
                        "<age>Ma_corrected_plus_dyntopo_diff_young_SW.nc files")
    p.add_argument("--ages", type=int, nargs="*", default=None,
                   help="explicit age list (default: all in --combined-dir)")
    p.add_argument("--fps", type=int, default=10)
    p.add_argument("--modes", nargs="+",
                   default=["combined", "dyntopo_diff"],
                   choices=["combined", "corrected", "dyntopo_diff",
                            "dyntopo_absolute"],
                   help="which video(s) to build: combined / corrected / "
                        "dyntopo_diff (per-step Δz) / dyntopo_absolute "
                        "(absolute past dyntopo on continents). "
                        "Default: combined dyntopo_diff.")
    p.add_argument("--keep-frames", action="store_true")
    p.add_argument("--force", action="store_true",
                   help="wipe cached per-frame PNGs under "
                        "video_frames_pygmt_SW_dyntopo/<mode>/ before rendering.")
    args = p.parse_args()

    if args.ages:
        ages = sorted(args.ages, reverse=True)
    else:
        ages = auto_discover_ages(args.combined_dir)
        if not ages:
            raise SystemExit(
                f"no combined NetCDFs in {args.combined_dir}. Run "
                f"`python build_dyntopo_diff_correction.py --source young "
                f"--ages $(seq 0 5 300) --step-myr 5` first.")
        ages = sorted(ages, reverse=True)
        print(f"auto-discovered {len(ages)} ages: "
              f"{ages[0]}..{ages[-1]} Ma")

    for mode in args.modes:
        sub = FRAME_DIR / mode
        sub.mkdir(exist_ok=True)
        existing = list(sub.glob("frame_*.png"))
        if args.force and existing:
            print(f"[{mode}] --force: wiping {len(existing)} cached frame(s)")
            for f in existing:
                f.unlink()
            existing = []
        elif existing and len(existing) > len(ages):
            print(f"[{mode}] {len(existing)} stale frames found "
                  f"(need {len(ages)}) — wiping folder")
            for f in existing:
                f.unlink()
            existing = []
        if existing and len(existing) >= len(ages):
            print(f"[{mode}] reusing {len(existing)} cached frames "
                  f"(pass --force to re-render)")
        for n, t in enumerate(ages):
            frame = sub / f"frame_{n:04d}.png"
            if frame.exists():
                continue
            print(f"[{mode}] rendering frame {n+1}/{len(ages)} at {t} Ma")
            render_frame(args.combined_dir, t, mode, frame)
        out_video = VIDEOS_DIR / f"SW_paleotopo_{mode}_{ages[0]}-{ages[-1]}Ma.mp4"
        if out_video.exists():
            out_video.unlink()
        stitch_video(sub / "frame_%04d.png", out_video,
                     args.fps, n_frames=len(ages))
        print(f"  wrote {out_video}")
        if not args.keep_frames:
            shutil.rmtree(sub)


if __name__ == "__main__":
    main()
