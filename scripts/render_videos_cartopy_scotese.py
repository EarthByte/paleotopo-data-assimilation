"""
=============================================================================
render_videos_cartopy_scotese.py  —  Preview videos for the S&W workflow
=============================================================================

Builds up to four MP4 videos from the Scotese & Wright corrected NetCDFs
using cartopy's Robinson projection (Winkel-Tripel-like).  This is the
*preview* renderer; for publication use render_videos_pygmt_scotese.py
which uses proper Winkel-Tripel via pyGMT.

OUTPUT VIDEOS  (paths_scotese.OUTPUT_DIR)
    SW_paleotopo_elevation_<oldest>-<youngest>Ma_preview.mp4
        ↳ corrected elevation
    SW_paleotopo_original_<oldest>-<youngest>Ma_preview.mp4
        ↳ original Scotese & Wright input
    SW_paleotopo_delta_<oldest>-<youngest>Ma_preview.mp4
        ↳ Δz = corrected − input
    SW_paleotopo_elevation_samples_<oldest>-<youngest>Ma_preview.mp4
        ↳ corrected + geochem control points (within ±Δt(t))

USAGE
    cd <project>/scripts
    python render_videos_cartopy_scotese.py                 # all four, 5 Myr
    python render_videos_cartopy_scotese.py --fps 12
    python render_videos_cartopy_scotese.py --modes elevation original
    python render_videos_cartopy_scotese.py --ages 0 50 100 200 500

OPTIONS
    --fps INT                 video frame rate (default 12)
    --ages INT [INT ...]      explicit age list (overrides auto)
    --cadence INT             render every N Ma (default: every S&W age)
    --modes ...               which video(s) to build
    --max-frames-per-call INT chunked-execution budget (for sandbox runs)

DEPENDENCIES  matplotlib, cartopy, netCDF4, ffmpeg

NOTES
    S&W has 107 distinct ages (5 Myr cadence with two gaps: 385 and 390 Ma
    are missing).  By default this script renders every available S&W age.
=============================================================================
"""
from __future__ import annotations
import argparse, subprocess, sys
from pathlib import Path
import numpy as np, netCDF4 as nc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths_scotese import CORRECTED_DIR, OUTPUT_DIR, CSV_PATH, VIDEOS_DIR
from sw_io import available_ages as sw_available_ages

OUT_DIR = OUTPUT_DIR
FRAME_DIR = OUTPUT_DIR / "video_frames_cartopy_SW"
FRAME_DIR.mkdir(parents=True, exist_ok=True)


# Adaptive temporal half-window — must match assimilate_scotese.dt_half()
def dt_half(t):
    if t < 200:  return 5.0
    if t < 500:  return 10.0
    if t < 800:  return 20.0
    return 30.0


def load(t):
    f = CORRECTED_DIR / f"{int(t)}Ma_corrected_SW.nc"
    with nc.Dataset(f) as d:
        lat = d.variables["lat"][:]; lon = d.variables["lon"][:]
        M = d.variables["M_orig"][:].astype(float)
        Mc = d.variables["M_corrected"][:].astype(float)
        delta = d.variables["delta"][:].astype(float)
        cont = d.variables["continent_mask"][:].astype(bool)
    return lat, lon, M, Mc, delta, cont


# ---------------------------------------------------------------------------
# Geochem cache + per-frame reconstruction (Scotese plate IDs via continental
# polygons; see sample_reconstruct_scotese.py).  Samples not in any
# continental polygon are dropped (would otherwise land in the ocean at
# older ages).
# ---------------------------------------------------------------------------
_GEOCHEM_CACHE = None
_SAMPLE_RECONSTRUCTOR = None
ELEV_MID_COLS = [
    "Isostatic_Elevation_absolute_km",
    "Brown_Isostatic_Elevation_absolute_km",
    "Davis_Isostatic_Elevation_absolute_km",
    "Condie_Isostatic_Elevation_absolute_km",
    "Herz_0.08_Isostatic_Elevation_absolute_km",
    "Herz_0.38 Isostatic_Elevation_absolute_km",
]
_PLATE_ID_CACHE = CSV_PATH.parent / "sample_plate_ids_SW.npy"


def _load_geochem():
    global _GEOCHEM_CACHE, _SAMPLE_RECONSTRUCTOR
    if _GEOCHEM_CACHE is not None:
        return _GEOCHEM_CACHE
    import pandas as pd
    from sample_reconstruct_scotese import ScoteseSampleReconstructor
    if _SAMPLE_RECONSTRUCTOR is None:
        _SAMPLE_RECONSTRUCTOR = ScoteseSampleReconstructor()

    df = pd.read_csv(CSV_PATH, low_memory=False)
    df = df.dropna(subset=["Lat", "Lon", "Age_Ma"])
    df = df[df.get("Missing_Risk_DeltaT (km)", 0).fillna(0) <= 5].copy()
    mid = df[ELEV_MID_COLS].astype(float)
    df["z_obs_m"] = mid.mean(axis=1, skipna=True) * 1000.0
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


def _samples_in_window(t, half=None, lat=None, lon=None, cont=None,
                       M=None, Mc=None):
    """Samples in ±half(t) window, reconstructed to age t via S&W rotations.

    If lat/lon/cont are provided, samples whose reconstructed position
    lands outside the per-slice continent mask are dropped — these come
    from terranes whose paleo-rotation has broken down at age t
    (polygon valid-time mismatch, missing terrane rotation, etc.).

    If M (original) and/or Mc (corrected) elevation arrays are also
    provided, samples landing on flooded cells (elevation ≤ 0 m) are
    additionally dropped.  This catches the rare case of a sample
    reconstructing into a drowned-interior basin or flooded shelf cell
    that is inside the polygon footprint but would visually plot over
    the ocean in the rendered map.  A sample is kept only if BOTH the
    original and corrected elevations at its cell are > 0 m (when both
    arrays are supplied)."""
    from sw_io import nearest_cell_index
    df = _load_geochem()
    if half is None: half = dt_half(t)
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


# Elevation range FIXED at −5000..+5000 m
# (so videos from both models can be compared side-by-side).
ELEVATION_RANGE = (-5000.0, 5000.0)
DATA_RANGE = {"elev": ELEVATION_RANGE, "delta": (-2000.0, 2000.0)}


def _nice_step(span):
    for s in (50, 100, 200, 250, 500, 1000):
        if span / s <= 50: return s
    return 2000.0


def compute_data_range(ages, modes):
    """Only computes the delta range from data — elevation is fixed."""
    needs_delta = "delta" in modes
    if not needs_delta:
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
    print(f"  delta cmap range derived from data: ±{d_abs:.0f} m  "
          f"(elevation range fixed at {ELEVATION_RANGE[0]:.0f} … {ELEVATION_RANGE[1]:.0f} m)")


def render_one(t, mode, frame):
    lat, lon, M, Mc, delta, cont = load(t)
    overlay_samples = False
    e_lo, e_hi = DATA_RANGE["elev"]; d_lo, d_hi = DATA_RANGE["delta"]
    elev_norm = mcolors.TwoSlopeNorm(vmin=e_lo, vcenter=0, vmax=e_hi)
    delta_norm = mcolors.TwoSlopeNorm(vmin=d_lo, vcenter=0, vmax=d_hi)
    if mode == "elevation":
        arr = Mc; cmap = plt.cm.terrain; norm = elev_norm
        title = "S&W corrected paleo-elevation"
        cb_label = "elevation (m)"
    elif mode == "elevation_samples":
        arr = Mc; cmap = plt.cm.terrain; norm = elev_norm
        title = (f"S&W corrected paleo-elevation + geochem control points "
                 f"(±{int(dt_half(t))} Myr)")
        cb_label = "elevation (m)"
        overlay_samples = True
    elif mode == "original":
        arr = M; cmap = plt.cm.terrain; norm = elev_norm
        title = "Scotese & Wright 2018 input paleo-elevation"
        cb_label = "elevation (m)"
    elif mode == "delta":
        # Mask Δz to the actual elevation-field paleo-land footprint
        # (M_orig ≥ 0 ∪ Mc ≥ 0), not the continent-polygon raster:
        # this removes the "ghost continent" Δz shadows in cells inside
        # the polygon footprint but below sea level in S&W18.
        land = (M >= 0) | (Mc >= 0)
        arr = np.where(land, delta, np.nan)
        cmap = plt.cm.RdBu_r; norm = delta_norm
        title = (f"Δz (corrected − S&W input) + geochem control points "
                 f"(±{int(dt_half(t))} Myr)")
        cb_label = "Δz (m)"
        overlay_samples = True
    else:
        raise ValueError(mode)

    fig = plt.figure(figsize=(11, 6))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson(central_longitude=0))
    ax.set_global()
    LON2D, LAT2D = np.meshgrid(lon, lat)
    pcm = ax.pcolormesh(LON2D, LAT2D, arr, cmap=cmap, norm=norm,
                        transform=ccrs.PlateCarree(), shading="auto",
                        rasterized=True)
    # No coastline contour on the videos.  At 1° / 5 Myr the contour
    # outlines pixel-edge artefacts that draw the eye away from the
    # actual elevation signal; the colour map itself defines the
    # coastline well enough at video resolution.  Static paper figures
    # keep the contour because they're read at a different scale.

    # No plate-boundary topology overlay on the S&W videos.  The S&W
    # plate model only resolves topologies for ages ≤ 100 Ma, so an
    # overlay would only appear in the last ~20 % of frames — not worth
    # the visual inconsistency.  The Merdith videos (which span the full
    # 0–1000 Ma with topologies) carry the topology overlay.

    if overlay_samples:
        samples = _samples_in_window(t, lat=lat, lon=lon, cont=cont,
                                     M=M, Mc=Mc)
        if len(samples):
            if mode == "delta":
                # Δz mode: small open black circles — position only.
                # Colour-coding samples by z_obs on a diverging Δz cmap
                # would be confusing (z_obs is absolute elevation, Δz is
                # a signed correction), so we just mark where the data
                # came from and let the underlying field show the
                # correction magnitude.
                ax.scatter(samples["rlon_t"].values, samples["rlat_t"].values,
                           facecolors="none", edgecolor="black",
                           s=18, linewidths=0.6,
                           transform=ccrs.PlateCarree(), zorder=5)
            else:
                # Elevation mode: filled circles coloured by z_obs on the
                # same terrain cmap as the background.
                ax.scatter(samples["rlon_t"].values, samples["rlat_t"].values,
                           c=samples["z_obs_m"].values,
                           cmap=cmap, norm=norm,
                           s=24, edgecolor="black", linewidths=0.6,
                           transform=ccrs.PlateCarree(), zorder=5)
        ax.text(0.02, 0.05, f"n = {len(samples)} samples",
                transform=ax.transAxes, fontsize=9,
                ha="left", va="bottom",
                bbox=dict(facecolor="white", alpha=0.7, edgecolor="none", pad=2))

    ax.gridlines(draw_labels=False, color="grey", lw=0.3, alpha=0.4)
    cax = fig.add_axes([0.18, 0.06, 0.64, 0.025])
    # extend='both' draws triangle arrows at both ends of the colorbar
    # showing that the data overflows the fixed ±5000 m range.
    cb = fig.colorbar(pcm, cax=cax, orientation="horizontal", extend="both")
    cb.set_label(cb_label, fontsize=10)
    ax.set_title(title, fontsize=14, weight="bold", pad=10)
    fig.text(0.95, 0.93, f"{int(t)} Ma", fontsize=18, weight="bold",
             ha="right", va="top", color="black")
    fig.subplots_adjust(left=0.02, right=0.98, top=0.94, bottom=0.13)
    fig.savefig(frame, dpi=140, bbox_inches="tight")
    plt.close(fig)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--cadence", type=int, default=None,
                   help="render every N Ma (default: all S&W ages)")
    p.add_argument("--fps", type=int, default=12)
    p.add_argument("--ages", type=int, nargs="*")
    p.add_argument("--modes", nargs="+",
                   default=["elevation", "delta", "original", "elevation_samples"],
                   choices=["elevation", "delta", "original", "elevation_samples"],
                   help="which video(s) to build (default: all four)")
    p.add_argument("--max-frames-per-call", type=int, default=0,
                   help="render at most N frames then exit (chunked runs)")
    p.add_argument("--force", action="store_true",
                   help="wipe cached per-frame PNGs (under video_frames_cartopy_SW/<mode>/) "
                        "and re-render from scratch. Use this after any code change to the "
                        "renderer or to the assimilation pipeline; otherwise existing frames "
                        "are reused and the stitched MP4 will silently contain stale visuals.")
    args = p.parse_args()

    if args.ages:
        ages = sorted(args.ages, reverse=True)
    else:
        all_ages = sw_available_ages()
        if args.cadence:
            ages = sorted([a for a in all_ages if a % args.cadence == 0], reverse=True)
        else:
            ages = sorted(all_ages, reverse=True)

    print(f"Scanning {len(ages)} NetCDFs to derive cmap range from actual data …")
    compute_data_range(ages, args.modes)
    rendered = 0
    for mode in args.modes:
        sub = FRAME_DIR / mode; sub.mkdir(exist_ok=True)
        existing = list(sub.glob("frame_*.png"))
        # --force wipes the per-mode folder unconditionally.
        if args.force and existing:
            print(f"[{mode}] --force: wiping {len(existing)} cached frame(s)")
            for f in existing:
                f.unlink()
            existing = []
        # Auto-wipe stale frames left over from earlier cadence/age-range runs.
        elif len(existing) > len(ages):
            print(f"[{mode}] {len(existing)} stale frames found (need {len(ages)}) — wiping folder")
            for f in existing:
                f.unlink()
            existing = []
        # Loud warning when frames are being reused — this hides upstream changes.
        if existing and len(existing) >= len(ages):
            print(f"[{mode}] reusing {len(existing)} cached frames "
                  f"(pass --force to re-render after a code change)")
        for n, t in enumerate(ages):
            frame = sub / f"frame_{n:04d}.png"
            if frame.exists(): continue
            if args.max_frames_per_call and rendered >= args.max_frames_per_call:
                print("hit per-call frame budget; exit clean for chunked driver")
                return
            render_one(t, mode, frame)
            rendered += 1
            if rendered % 20 == 0:
                print(f"[{mode}] rendered {rendered} frames so far (latest age {t} Ma)")
        existing = sorted(sub.glob("frame_*.png"))
        if len(existing) == len(ages):
            out = VIDEOS_DIR / f"SW_paleotopo_{mode}_{ages[0]}-{ages[-1]}Ma_preview.mp4"
            if out.exists(): out.unlink()
            cmd = ["ffmpeg", "-y", "-framerate", str(args.fps),
                   "-i", str(sub / "frame_%04d.png"),
                   "-frames:v", str(len(ages)),
                   "-c:v", "libx264", "-pix_fmt", "yuv420p",
                   "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2",
                   "-crf", "20", str(out)]
            print(" ".join(cmd))
            subprocess.run(cmd, check=True)
            print(f"  wrote {out}")


if __name__ == "__main__":
    main()
