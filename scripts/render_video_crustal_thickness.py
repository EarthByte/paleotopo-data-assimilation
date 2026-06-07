"""
=============================================================================
render_video_crustal_thickness.py
=============================================================================

Stitches a Robinson-projection MP4 sweeping through every Phanerozoic
slice (0–540 Ma) of the Airy-isostasy crustal-thickness field z_c
produced by `derive_crustal_thickness.py`.

INPUTS  (must exist; produced by derive_crustal_thickness.py --all)
    data/corrected/<age>Ma_crustal_thickness_SW.nc

OUTPUTS
    outputs/video_frames_crustal_thickness_SW/frame_NNNN.png   (cache)
    outputs/SW_crustal_thickness_540-0Ma.mp4                   (final)

USAGE
    cd <repo root>/scripts
    python render_video_crustal_thickness.py                 # all ages
    python render_video_crustal_thickness.py --fps 8
    python render_video_crustal_thickness.py --cadence 10    # every 10 Ma
    python render_video_crustal_thickness.py --force         # re-render frames
=============================================================================
"""
from __future__ import annotations
import argparse, subprocess, sys
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
from paths_scotese import CORRECTED_DIR, OUTPUT_DIR

FRAME_DIR = OUTPUT_DIR / "video_frames_crustal_thickness_SW"
OUT_DIR   = OUTPUT_DIR
FRAME_DIR.mkdir(parents=True, exist_ok=True)
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Crustal-thickness colour scale.  20 km lower bound captures rifted
# margins; 70 km upper bound captures collision-belt cores.  Values
# outside the range render as extended colorbar triangles.
Z_LO, Z_HI = 20.0, 70.0


def _make_haxby_cmap():
    """Standard GMT 'haxby' palette as a matplotlib colormap.

    Haxby is a perceptually-uniform purple→blue→green→yellow→orange→
    red→brown ramp originally designed for ocean bathymetry but well
    suited for any unsigned scalar field.  Stops are the 11 canonical
    GMT haxby RGB nodes; matplotlib linearly interpolates between them.
    """
    stops = np.array([
        (  9,  60, 168),  # deep purple
        ( 34, 100, 196),
        ( 58, 142, 226),
        (106, 191, 234),
        (153, 224, 230),
        (189, 222, 192),
        (213, 213, 132),
        (239, 209,  61),
        (237, 167,  43),
        (223, 109,  27),
        (200,  56,  13),  # dark red-brown
    ], dtype=float) / 255.0
    return mcolors.LinearSegmentedColormap.from_list("gmt_haxby", stops, N=512)


CMAP = _make_haxby_cmap()
NORM = mcolors.Normalize(vmin=Z_LO, vmax=Z_HI)


def _available_ages():
    files = sorted(CORRECTED_DIR.glob("*Ma_crustal_thickness_SW.nc"),
                   key=lambda f: int(f.name.split("Ma")[0]))
    return [int(f.name.split("Ma")[0]) for f in files]


def render_frame(age_ma: int, out_path: Path):
    f = CORRECTED_DIR / f"{int(age_ma)}Ma_crustal_thickness_SW.nc"
    with nc.Dataset(f) as d:
        lat = np.asarray(d.variables["lat"][:])
        lon = np.asarray(d.variables["lon"][:])
        z_c = np.asarray(d.variables["z_c"][:], dtype=float)
        # h_m is the M_corrected echo stored alongside z_c — we use it
        # to draw the sea-level coastline contour so the outline matches
        # the elevation grid the thickness was derived from, rather
        # than the continent-polygon raster (which is offset before
        # ~350 Ma due to a polygon/rotation mismatch).
        h_m = np.asarray(d.variables["h_m"][:], dtype=float)

    fig = plt.figure(figsize=(11, 6))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson(central_longitude=0))
    ax.set_global()
    LON2D, LAT2D = np.meshgrid(lon, lat)
    pcm = ax.pcolormesh(LON2D, LAT2D, z_c, cmap=CMAP, norm=NORM,
                        transform=ccrs.PlateCarree(),
                        shading="auto", rasterized=True)
    # No coastline contour — at 1° / 5 Myr it traces pixel-edge
    # artefacts more than the actual paleo-coastline.
    ax.gridlines(lw=0.25, color="grey", alpha=0.4)
    ax.set_title(
        "Derived continental crustal thickness "
        "(Airy isostasy on corrected S&W18 elevation)",
        fontsize=13, weight="bold", pad=10
    )
    fig.text(0.95, 0.93, f"{int(age_ma)} Ma", fontsize=18, weight="bold",
             ha="right", va="top", color="black")

    cax = fig.add_axes([0.18, 0.06, 0.64, 0.025])
    cb = fig.colorbar(pcm, cax=cax, orientation="horizontal", extend="both")
    cb.set_label("z_c (km)", fontsize=10)
    fig.subplots_adjust(left=0.02, right=0.98, top=0.94, bottom=0.13)
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def stitch(n_frames: int, fps: int, age_lo: int, age_hi: int):
    out = OUT_DIR / f"SW_crustal_thickness_{age_hi}-{age_lo}Ma.mp4"
    if out.exists():
        out.unlink()
    cmd = ["ffmpeg", "-y", "-framerate", str(fps),
           "-i", str(FRAME_DIR / "frame_%04d.png"),
           "-frames:v", str(n_frames),
           "-c:v", "libx264", "-pix_fmt", "yuv420p",
           "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2",
           "-crf", "20", str(out)]
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)
    print(f"  wrote {out}")
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--fps", type=int, default=6,
                   help="frames per second (default 6 — ~18 s for the full "
                        "Phanerozoic at the native 5 Myr cadence)")
    p.add_argument("--cadence", type=int, default=None,
                   help="render every N Ma (default: every available age)")
    p.add_argument("--ages", type=int, nargs="*",
                   help="explicit list of ages (overrides --cadence)")
    p.add_argument("--force", action="store_true",
                   help="wipe cached frames before rendering")
    args = p.parse_args()

    if args.ages:
        ages = sorted(args.ages, reverse=True)
    else:
        all_ages = _available_ages()
        if not all_ages:
            sys.exit("no crustal-thickness NetCDFs found in "
                     f"{CORRECTED_DIR} — run derive_crustal_thickness.py --all first")
        if args.cadence:
            ages = sorted([a for a in all_ages if a % args.cadence == 0],
                          reverse=True)
        else:
            ages = sorted(all_ages, reverse=True)

    existing = list(FRAME_DIR.glob("frame_*.png"))
    # The cache is indexed by frame number, not by age — so frame_0000.png
    # produced by an earlier run at a different cadence is at a DIFFERENT
    # age than frame_0000 of this run.  The only safe reuse is when the
    # existing-frame count matches the new age list EXACTLY (same cadence,
    # same age range).  Any mismatch wipes.
    if args.force and existing:
        print(f"--force: wiping {len(existing)} cached frame(s)")
        for f in existing:
            f.unlink()
        existing = []
    elif existing and len(existing) != len(ages):
        print(f"{len(existing)} cached frames don't match new run "
              f"(need {len(ages)}) — wiping")
        for f in existing:
            f.unlink()
        existing = []
    elif existing and len(existing) == len(ages):
        print(f"reusing {len(existing)} cached frames "
              f"(pass --force to re-render after a code change)")

    print(f"rendering {len(ages)} frames at {args.fps} fps "
          f"(ages {ages[0]} → {ages[-1]} Ma) …")
    for n, t in enumerate(ages):
        frame = FRAME_DIR / f"frame_{n:04d}.png"
        if frame.exists():
            continue
        render_frame(t, frame)
        if (n + 1) % 20 == 0:
            print(f"  rendered {n + 1}/{len(ages)} frames (latest age {t} Ma)")

    stitch(len(ages), args.fps, age_lo=ages[-1], age_hi=ages[0])


if __name__ == "__main__":
    main()
