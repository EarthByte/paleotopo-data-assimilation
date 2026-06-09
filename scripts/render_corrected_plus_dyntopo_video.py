"""
=============================================================================
render_corrected_plus_dyntopo_video.py  —  Preview videos of corrected-SW +
                                           Young 2022 dynamic topography
=============================================================================

Sibling of render_videos_cartopy_scotese.py. Reads the combined NetCDFs
produced by add_dyntopo_to_corrected_scotese.py and renders Robinson-
projection MP4 previews.

OUTPUT VIDEOS  (under OUTPUT_DIR)
    SW_paleotopo_combined_<oldest>-<youngest>Ma_preview.mp4
        ↳ corrected paleo-elevation + Young 2022 dyntopo
          (M_combined), terrain colourmap (TwoSlopeNorm at 0)

    SW_paleotopo_dyntopo_<oldest>-<youngest>Ma_preview.mp4
        ↳ Young 2022 dyntopo time-difference (z_dyntopo_diff) only,
          RdBu_r at ±1500 m

    SW_paleotopo_corrected_dt_<oldest>-<youngest>Ma_preview.mp4
        ↳ corrected paleo-elevation (M_corrected, for reference),
          same age window as the dyntopo videos

Only the [--age-min, --age-max] window covered by both corrected-SW and
the rotated dyntopo grids is rendered (0..250 Ma by default; that's the
Scotese ↔ Merdith plate-motion convergence window).

USAGE
    cd <project>/scripts_Scotese
    python add_dyntopo_to_corrected_scotese.py --dyntopo-dir <dyntopo>
    python render_corrected_plus_dyntopo_video.py                    # all three videos
    python render_corrected_plus_dyntopo_video.py --modes combined dyntopo
    python render_corrected_plus_dyntopo_video.py --fps 6 --age-max 250

OPTIONS
    --combined-dir PATH       directory of <age>Ma_corrected_plus_dyntopo_SW.nc
                              (default: PROJECT_ROOT/data/corrected_Scotese_plus_dyntopo)
    --age-min INT             default 0
    --age-max INT             default 250 (caps to dyntopo availability)
    --cadence-myr INT         default 5
    --fps INT                 default 12
    --modes ...               which video(s) to build
    --vlim-dyntopo FLOAT      ±m for the dyntopo cmap (default 1500)
    --force                   re-render all cached PNG frames

DEPENDENCIES  matplotlib, cartopy, netCDF4, ffmpeg
=============================================================================
"""
from __future__ import annotations
import argparse, re, subprocess, sys
from pathlib import Path

import numpy as np
import netCDF4 as nc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths_scotese import OUTPUT_DIR, PROJECT_ROOT

COMBINED_FNAME_RE = re.compile(r"^(\d+)Ma_corrected_plus_dyntopo_diff_young_SW\.nc$")
COMBINED_FNAME_FMT = "{age}Ma_corrected_plus_dyntopo_diff_young_SW.nc"

OUT_DIR = OUTPUT_DIR
FRAME_DIR = OUTPUT_DIR / "video_frames_cartopy_SW_dyntopo"
FRAME_DIR.mkdir(parents=True, exist_ok=True)

# Match the elevation range used in render_videos_cartopy_scotese.py so
# the combined-topo video can be compared side-by-side with the corrected
# video at the same vmin/vmax.
ELEVATION_RANGE = (-5000.0, 5000.0)
DEFAULT_DYNTOPO_VLIM = 1500.0


def _load(combined_dir: Path, age: int) -> tuple:
    """Read M_corrected, M_combined, z_dyntopo_diff and the continent mask
    for `age` from the build_dyntopo_diff_correction.py output NetCDF.

    The new diff-NetCDF format stores the plate-frame-computed +
    Scotese-paleomag-rotated dyntopo difference as `z_dyntopo_diff` (the
    correction signal actually added to M_corrected) rather than the
    absolute past dyntopo, and does NOT carry a continent mask — we
    derive it from the corresponding corrected-S&W file instead."""
    fp = combined_dir / COMBINED_FNAME_FMT.format(age=age)
    with nc.Dataset(fp) as d:
        lat = np.asarray(d["lat"][:], dtype=np.float64)
        lon = np.asarray(d["lon"][:], dtype=np.float64)
        Mc  = np.asarray(d["M_corrected"][:], dtype=float)
        z   = np.asarray(d["z_dyntopo_diff"][:], dtype=float)
        Mcomb = np.asarray(d["M_combined"][:], dtype=float)
    corrected_sw = (PROJECT_ROOT / "data" / "corrected_Scotese"
                    / f"{int(age)}Ma_corrected_SW.nc")
    if corrected_sw.exists():
        with nc.Dataset(corrected_sw) as d:
            cont = np.asarray(d["continent_mask"][:], dtype=bool)
    else:
        cont = (Mc >= 0)
    return lat, lon, Mc, z, Mcomb, cont


def _available_ages(combined_dir: Path) -> list[int]:
    if not combined_dir.is_dir():
        return []
    return sorted(
        int(COMBINED_FNAME_RE.match(p.name).group(1))
        for p in combined_dir.iterdir()
        if COMBINED_FNAME_RE.match(p.name)
    )


def render_one(t: int, mode: str, frame_path: Path, combined_dir: Path,
               vlim_dyntopo: float) -> None:
    lat, lon, Mc, z, Mcomb, cont = _load(combined_dir, t)
    e_lo, e_hi = ELEVATION_RANGE
    elev_norm = mcolors.TwoSlopeNorm(vmin=e_lo, vcenter=0.0, vmax=e_hi)
    dyn_norm  = mcolors.TwoSlopeNorm(vmin=-vlim_dyntopo, vcenter=0.0, vmax=vlim_dyntopo)

    if mode == "combined":
        arr = Mcomb; cmap = plt.cm.terrain; norm = elev_norm
        title = "S&W corrected paleo-elevation + Young 2022 dyntopo"
        cb_label = "elevation (m)"
    elif mode == "corrected":
        arr = Mc; cmap = plt.cm.terrain; norm = elev_norm
        title = "S&W corrected paleo-elevation"
        cb_label = "elevation (m)"
    elif mode == "dyntopo":
        arr = z; cmap = plt.cm.RdBu_r; norm = dyn_norm
        title = "Young 2022 dynamic topography (Scotese frame)"
        cb_label = "dyntopo (m)"
    else:
        raise ValueError(mode)

    fig = plt.figure(figsize=(11, 6))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson(central_longitude=0))
    ax.set_global()
    LON2D, LAT2D = np.meshgrid(lon, lat)
    pcm = ax.pcolormesh(LON2D, LAT2D, arr, cmap=cmap, norm=norm,
                        transform=ccrs.PlateCarree(), shading="auto",
                        rasterized=True)
    ax.contour(LON2D, LAT2D, cont.astype(float),
               levels=[0.5], colors="black", linewidths=0.4,
               transform=ccrs.PlateCarree())
    ax.gridlines(draw_labels=False, color="grey", lw=0.3, alpha=0.4)
    cax = fig.add_axes([0.18, 0.06, 0.64, 0.025])
    cb = fig.colorbar(pcm, cax=cax, orientation="horizontal", extend="both")
    cb.set_label(cb_label, fontsize=10)
    ax.set_title(title, fontsize=14, weight="bold", pad=10)
    fig.text(0.95, 0.93, f"{int(t)} Ma", fontsize=18, weight="bold",
             ha="right", va="top", color="black")
    fig.subplots_adjust(left=0.02, right=0.98, top=0.94, bottom=0.13)
    fig.savefig(frame_path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--combined-dir", type=Path,
                   default=PROJECT_ROOT / "data" / "corrected_Scotese_plus_dyntopo",
                   help="directory of <age>Ma_corrected_plus_dyntopo_SW.nc files")
    p.add_argument("--age-min", type=int, default=0)
    p.add_argument("--age-max", type=int, default=250)
    p.add_argument("--cadence-myr", type=int, default=5)
    p.add_argument("--fps", type=int, default=12)
    p.add_argument("--vlim-dyntopo", type=float, default=DEFAULT_DYNTOPO_VLIM,
                   help="symmetric colour range for the dyntopo cmap; default ±1500 m")
    p.add_argument("--modes", nargs="+",
                   default=["combined", "dyntopo", "corrected"],
                   choices=["combined", "dyntopo", "corrected"],
                   help="which video(s) to build (default: all three)")
    p.add_argument("--force", action="store_true",
                   help="wipe cached per-frame PNGs (under FRAME_DIR/<mode>/) "
                        "and re-render. Use after any code change to the "
                        "renderer; otherwise existing frames are reused.")
    p.add_argument("--max-frames-per-call", type=int, default=0,
                   help="render at most N frames then exit (chunked runs)")
    args = p.parse_args()

    avail = _available_ages(args.combined_dir)
    if not avail:
        raise SystemExit(
            f"no <age>Ma_corrected_plus_dyntopo_SW.nc files under "
            f"{args.combined_dir}. Run add_dyntopo_to_corrected_scotese.py "
            f"first."
        )
    requested = list(range(args.age_min,
                           args.age_max + 1,
                           args.cadence_myr))
    ages = [a for a in requested if a in avail]
    missing = [a for a in requested if a not in avail]
    if missing:
        print(f"WARNING: skipping {len(missing)} requested ages not present in "
              f"combined-dir: {missing[:10]}{'...' if len(missing) > 10 else ''}")
    ages = sorted(ages, reverse=True)  # render oldest → youngest, like the SW video
    if not ages:
        raise SystemExit("no ages to render after intersection with available files")
    print(f"will render {len(ages)} frames per mode "
          f"({ages[0]} → {ages[-1]} Ma at {args.cadence_myr} Myr cadence)")

    rendered = 0
    for mode in args.modes:
        sub = FRAME_DIR / mode; sub.mkdir(exist_ok=True)
        existing = list(sub.glob("frame_*.png"))
        if args.force and existing:
            print(f"[{mode}] --force: wiping {len(existing)} cached frame(s)")
            for f in existing:
                f.unlink()
            existing = []
        elif len(existing) > len(ages):
            print(f"[{mode}] {len(existing)} stale frames found (need {len(ages)}) — wiping folder")
            for f in existing:
                f.unlink()
            existing = []
        if existing and len(existing) >= len(ages):
            print(f"[{mode}] reusing {len(existing)} cached frames "
                  f"(pass --force to re-render after a code change)")
        for n, t in enumerate(ages):
            frame = sub / f"frame_{n:04d}.png"
            if frame.exists():
                continue
            if args.max_frames_per_call and rendered >= args.max_frames_per_call:
                print("hit per-call frame budget; exit clean for chunked driver")
                return 0
            render_one(t, mode, frame, args.combined_dir, args.vlim_dyntopo)
            rendered += 1
            if rendered % 20 == 0:
                print(f"[{mode}] rendered {rendered} frames so far (latest age {t} Ma)")
        existing = sorted(sub.glob("frame_*.png"))
        if len(existing) == len(ages):
            tag = {"combined": "combined",
                   "corrected": "corrected_dt",
                   "dyntopo": "dyntopo"}[mode]
            out = OUT_DIR / f"SW_paleotopo_{tag}_{ages[0]}-{ages[-1]}Ma_preview.mp4"
            if out.exists():
                out.unlink()
            cmd = ["ffmpeg", "-y", "-framerate", str(args.fps),
                   "-i", str(sub / "frame_%04d.png"),
                   "-frames:v", str(len(ages)),
                   "-c:v", "libx264", "-pix_fmt", "yuv420p",
                   "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2",
                   "-crf", "20", str(out)]
            print(" ".join(cmd))
            subprocess.run(cmd, check=True)
            print(f"  wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
