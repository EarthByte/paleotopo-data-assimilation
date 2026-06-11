"""
=============================================================================
build_dyntopo_diff_correction.py  —  PER-STEP dyntopo increment correction
                                      (plate-frame, paleomag-rotated)
=============================================================================

Builds <T>Ma_corrected_plus_dyntopo_diff_<source>_SW.nc files by adding a
PER-STEP increment of plate-frame dynamic topography (Young 2022 OR
Braz 2021) to the geochem-corrected Scotese & Wright paleo-DEM:

    Δz(t)        = dyntopo(t) − dyntopo(t − Δt)   (Δt = --step-myr, default 5)
    M_combined(t) = M_corrected(t) + Δz_paleomag(t)

Each age t therefore inherits the dyntopo CONTRIBUTION to topographic change
over the most recent Δt-Myr interval, NOT a cumulative-from-present
correction.  The 0 Ma slice serves as a reference for the t = Δt step only;
at every later step the reference SLIDES one step back (t = 10 Ma uses
5 Ma; t = 15 uses 10; etc.).  At t = 0 the correction is identically zero
by construction, so M_combined(0) = M_corrected(0) = observed topography.

The intent is to visualise where mantle flow is most actively reshaping the
surface at each frame — not to reconstruct what cumulative dyntopo applied
on top of M_corrected would look like.

WHY PLATE FRAME?  In mantle reference frame, continents move over time.
Subtracting two mantle-frame dyntopo grids at the same lat/lon cell compares
the dyntopo under different parts of different continents at different times
— physically meaningless. In plate frame, every grid cell is rigidly attached
to its continent, so dyntopo(t) − dyntopo(t − Δt) gives the time-change at
that specific point on that specific continent over the Δt interval.

SIGN CONVENTION
  dyntopo > 0  → surface is uplifted by mantle dynamics at that cell
  Δz(t) > 0    → that cell was uplifted (or experienced less subsidence)
                  by mantle dynamics in the Δt interval ending at age t
  ⇒ M_combined(t) = M_corrected(t)
                    + [dyntopo(t) − dyntopo(t − Δt)] rotated to paleomag(t)

INPUTS
  Young 2022 plate-frame grids:
    data/Young2022_gld428_grids_20Myr/gld428_PlateFrameGrid<age>.nc
    20-Myr cadence (0, 20, 40, ..., 980).  Variable: 'z' (or similar — sniffed).

  Braz 2021 gmcm9 plate-frame grids:
    data/Braz_etal_2021-_dynatopo_and_rotations/gmcm9_plate_ref_frame_grids/
    <age>.00.nc  for 0, 4, 9, 14, 19, ..., 99 (5-Myr cadence, offset 1 Myr)
    plus 104, 109, ..., 150 Ma.

  Geochem-corrected paleo-DEM (input to add the dyntopo diff to):
    Paleotopo_data_assimilation/data/corrected_Scotese/<age>Ma_corrected_SW.nc

  Scotese plate model (for partition + rotation):
    Scotese_Wright_PlateModel_2023.rot
    Scotese_Wright_PresentDay_ContinentalPolygons.gpml

OUTPUTS
  Paleotopo_data_assimilation/data/corrected_Scotese_plus_dyntopo_diff_<source>/
    <T>Ma_corrected_plus_dyntopo_diff_<source>_SW.nc
  with M_corrected (input), z_dyntopo_diff (rotated to Scotese frame), and
  M_combined (sum, land>0 guarded) — all CF-1.10 + ACDD-1.3 compliant.

DEPENDENCIES
  gplately ≥ 1.0  (already in envi_gospl per task #2)
  pygplates ≥ 1.0
  netCDF4, xarray, numpy
  netcdf_io.py — CF writer helper (sibling file in this directory)

USAGE
  python scripts_Scotese/build_dyntopo_diff_correction.py \\
      --source young \\
      --ages 5 15 35 55
=============================================================================
"""
from __future__ import annotations
import argparse
import os
import re
import sys
from pathlib import Path

import numpy as np
import xarray as xr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from netcdf_io import write_cf_grid, git_describe

# Resolve paleotopo root from this script's location (scripts_Scotese/.. → repo).
PALEOTOPO_ROOT = HERE.parent

# Where the Young 2022 / Braz 2021 plate-frame dyntopo grids live.
# Default points at the project's own data/ directory; override via
# env var PALEOTOPO_DATA_ROOT or the --young-dir / --braz-dir CLI flags.
DATA_ROOT = Path(os.environ.get(
    "PALEOTOPO_DATA_ROOT",
    str(PALEOTOPO_ROOT / "data")))


# ---------------------------------------------------------------------------
# Source-specific configuration
# ---------------------------------------------------------------------------

YOUNG_DIR = DATA_ROOT / "Young2022_gld428_grids_20Myr"
BRAZ_DIR  = (DATA_ROOT / "Braz_etal_2021-_dynatopo_and_rotations"
             / "gmcm9_plate_ref_frame_grids")

CORRECTED_DIR = PALEOTOPO_ROOT / "data" / "corrected_Scotese"

# Scotese 2023 plate model for partition + rotation
SCOTESE_PLATE_DIR = (PALEOTOPO_ROOT / "data" /
                     "Scotese_Wright_plate_model_revised")
SCOTESE_ROT_FILE = SCOTESE_PLATE_DIR / "Scotese_Wright_PlateModel_2023.rot"
SCOTESE_POLYGONS = SCOTESE_PLATE_DIR / "Scotese_Wright_PresentDay_ContinentalPolygons.gpml"

# Sea-level floor for the dyntopo-correction guard. When dyntopo(T) − dyntopo(0)
# pushes a cell that was originally above sea level (M_orig > 0) below 0 m
# in M_combined, we clip it to LAND_MIN_ELEV_M. Setting this to 0.0 is the
# "sea level" interpretation — the region is now exactly at the shoreline,
# physically reasonable for genuine dyntopo-driven coastal inundation.
# Higher values would artificially keep the cell as "definitely land" which
# is dishonest about what dyntopo did.
LAND_MIN_ELEV_M = 0.0


# ---------------------------------------------------------------------------
# Source-grid resolution + loading
# ---------------------------------------------------------------------------

def _young_age_to_file(age_ma: float) -> Path:
    """Young 2022 file for the nearest 20-Myr keyframe, or two for interpolation."""
    return YOUNG_DIR / f"gld428_PlateFrameGrid{int(round(age_ma))}.nc"


def _braz_age_to_file(age_ma: float) -> Path:
    """Braz 2021 file: cadence is 0, 4, 9, 14, 19, 24, 29, 34, 39, 44, 49, 54, 59, 64, 69, ...
    Snap to nearest available age."""
    available = sorted([
        int(p.stem.split(".")[0]) for p in BRAZ_DIR.iterdir()
        if re.match(r"^\d+\.00\.nc$", p.name)
    ])
    snapped = min(available, key=lambda a: abs(a - age_ma))
    return BRAZ_DIR / f"{snapped:d}.00.nc"


def load_dyntopo_at_age(source: str, age_ma: float) -> tuple[xr.DataArray, int]:
    """Return (dyntopo DataArray on plate-frame lat/lon, actual_source_age).
    For Young: linearly interpolate between the two flanking 20-Myr keyframes
    if age_ma isn't exactly on a keyframe. For Braz: snap to nearest 5-Myr-offset.
    """
    if source == "young":
        # Linear interpolation between the two flanking 20-Myr keyframes
        lo = int(np.floor(age_ma / 20) * 20)
        hi = lo + 20
        fp_lo = _young_age_to_file(lo); fp_hi = _young_age_to_file(hi)
        if not fp_lo.exists():
            raise FileNotFoundError(fp_lo)
        if not fp_hi.exists():
            raise FileNotFoundError(fp_hi)
        ds_lo = xr.open_dataset(fp_lo, decode_times=False)
        ds_hi = xr.open_dataset(fp_hi, decode_times=False)
        var = _sniff_dyntopo_var(ds_lo)
        w = (age_ma - lo) / 20.0   # 0 → lo, 1 → hi
        da = (1 - w) * ds_lo[var] + w * ds_hi[var]
        # Carry over coordinates from one of them
        da = da.assign_coords(ds_lo[var].coords)
        return _normalise_lat_lon(da), int(round(age_ma))
    elif source == "braz":
        fp = _braz_age_to_file(age_ma)
        ds = xr.open_dataset(fp, decode_times=False)
        var = _sniff_dyntopo_var(ds)
        actual_age = int(fp.stem.split(".")[0])
        return _normalise_lat_lon(ds[var]), actual_age
    else:
        raise ValueError(f"unknown source: {source}")


def _sniff_dyntopo_var(ds: xr.Dataset) -> str:
    """Find the dyntopo data variable in a dyntopo NetCDF. Common names: z,
    dt, topo, dyntopo. Falls back to the first non-coordinate float variable."""
    for c in ("z", "dt", "topo", "dyntopo", "elevation"):
        if c in ds.data_vars:
            return c
    for v, da in ds.data_vars.items():
        if np.issubdtype(da.dtype, np.floating) and da.ndim == 2:
            return v
    raise KeyError(f"No 2-D float variable in {ds.data_vars}")


def _normalise_lat_lon(da: xr.DataArray) -> xr.DataArray:
    """Rename y/x→lat/lon, force ascending lat, lon in [-180, +180]."""
    rename = {}
    for src, tgt in (("latitude", "lat"), ("longitude", "lon"),
                     ("y", "lat"), ("x", "lon")):
        if src in da.coords:
            rename[src] = tgt
    if rename:
        da = da.rename(rename)
    if da["lat"].values[0] > da["lat"].values[-1]:
        da = da.isel(lat=slice(None, None, -1))
    if float(da["lon"].max()) > 180:
        da = da.assign_coords(lon=(((da["lon"] + 180) % 360) - 180))
        da = da.sortby("lon")
    return da


def load_corrected_scotese(age_ma: int) -> tuple[xr.DataArray, xr.DataArray]:
    """Load (M_corrected, M_orig) at the target age from the geochem-corrected
    Scotese directory."""
    fp = CORRECTED_DIR / f"{age_ma}Ma_corrected_SW.nc"
    if not fp.exists():
        raise FileNotFoundError(fp)
    ds = xr.open_dataset(fp, decode_times=False)
    return _normalise_lat_lon(ds["M_corrected"]), _normalise_lat_lon(ds["M_orig"])


# ---------------------------------------------------------------------------
# Plate-frame → Scotese paleomag frame via gplately
# ---------------------------------------------------------------------------

def cookie_cut_and_rotate(plate_frame: xr.DataArray, age_ma: float,
                          rot_file: Path, polygons_file: Path) -> xr.DataArray:
    """Take a present-day-position plate-frame raster (continents in their
    present-day positions) and partition + rotate each continental block to
    its position at `age_ma` using the Scotese 2023 plate model. Returns a
    raster on the same lat/lon grid but with values rotated to paleomag frame.

    Uses gplately.Raster.reconstruct(reverse=True) — "reverse" because we
    want present-day → past, not past → present-day. (Standard gplately
    reconstruction is from a past raster TO present-day; we want the inverse.)
    """
    import gplately
    import pygplates
    # gplately ≥ 1.0 PlateReconstruction signature: rotation_model + optional
    # topology_features + static_polygons. The rotation_model can be a path,
    # a FeatureCollection, or a pygplates.RotationModel.
    plate_model = gplately.PlateReconstruction(
        rotation_model=pygplates.RotationModel(str(rot_file)),
        static_polygons=pygplates.FeatureCollection(str(polygons_file)),
    )
    lat = plate_frame["lat"].values
    lon = plate_frame["lon"].values
    raster = gplately.Raster(
        data=plate_frame.values.astype(np.float32),
        plate_reconstruction=plate_model,
        time=0,    # present-day
        extent=[float(lon.min()), float(lon.max()),
                float(lat.min()), float(lat.max())],
    )
    # reconstruct to past age — equivalent to "where would today's continents
    # have been at age_ma". The static_polygons attached to plate_model are
    # used for partitioning; no partitioning_features kwarg needed.
    reconstructed = raster.reconstruct(
        time=float(age_ma),
        threads=1,
        anchor_plate_id=0,
        inplace=False,
    )
    # gplately.Raster.reconstruct returns a new Raster (or None if inplace).
    # Resolve to the underlying numpy array on the same lat/lon target grid.
    if reconstructed is None:
        data_out = np.asarray(raster.data)
    else:
        data_out = np.asarray(reconstructed.data)
    out = xr.DataArray(
        data=data_out,
        dims=("lat", "lon"),
        coords={"lat": lat, "lon": lon},
        name="z_dyntopo_diff_scotese",
    )
    return out


# ---------------------------------------------------------------------------
# Main per-age processor
# ---------------------------------------------------------------------------

def process_age(target_age_ma: int, source: str, *,
                rot_file: Path, polygons_file: Path,
                out_dir: Path, step_myr: int = 5,
                verbose: bool = True) -> Path:
    """At age t, the dyntopo correction is the PER-STEP increment

        Δz(t) = z_dyntopo(t) − z_dyntopo(t − Δt),    Δt = `step_myr`

    rather than the cumulative-from-present z(t) − z(0).  At t = 0 the
    correction is identically zero by construction (no predecessor age),
    so M_combined(0) = M_corrected(0) = observed topography.  The
    reference age slides with t so that 0 Ma serves as a reference for
    the t = Δt step only, never thereafter.
    """
    if verbose:
        print(f"\n=== {target_age_ma} Ma — dyntopo source: {source}, "
              f"step Δt = {step_myr} Myr ===")

    ref_age_ma = max(0, target_age_ma - step_myr)

    if target_age_ma == 0:
        # No predecessor — Δz is identically zero, M_combined ≡ M_corrected.
        # Load M_corrected just so the downstream cookie-cut/rotate step has
        # a target grid to project zeros onto.
        M_corr_pre, _ = load_corrected_scotese(target_age_ma)
        diff_plate = xr.zeros_like(M_corr_pre)
        diff_plate.name = "z"
        src_age_T = src_age_ref = 0
        if verbose:
            print("  t = 0: Δz set to zero by construction")
    else:
        # 1. dyntopo at t in plate frame
        dyntopo_T, src_age_T = load_dyntopo_at_age(source, target_age_ma)
        if verbose:
            print(f"  dyntopo({target_age_ma} Ma)  ← {source} source age "
                  f"{src_age_T} Ma")
        # 2. dyntopo at t − Δt in plate frame (per-step reference, NOT 0 Ma)
        dyntopo_ref, src_age_ref = load_dyntopo_at_age(source, float(ref_age_ma))
        if verbose:
            print(f"  dyntopo({ref_age_ma} Ma)  ← {source} source age "
                  f"{src_age_ref} Ma  (Δt={step_myr} Myr reference)")
        # 3. PER-STEP difference in plate frame: Δz = z(t) − z(t − Δt)
        diff_plate = dyntopo_T - dyntopo_ref
        if verbose:
            print(f"  Δz_plate range: [{float(diff_plate.min()):+.0f}, "
                  f"{float(diff_plate.max()):+.0f}] m")

    # 4-6. cookie-cut + rotate to Scotese paleomag frame
    if verbose:
        print(f"  rotating into Scotese paleomag frame ...")
    diff_scotese = cookie_cut_and_rotate(diff_plate, target_age_ma,
                                          rot_file, polygons_file)
    if verbose:
        print(f"  diff_scotese range: [{float(diff_scotese.min()):+.0f}, "
              f"{float(diff_scotese.max()):+.0f}] m")

    # 7. load M_corrected at target age + add the diff
    M_corr, M_orig = load_corrected_scotese(target_age_ma)
    # Make sure all three are on the same grid; interpolate diff if needed
    if not (np.array_equal(diff_scotese["lat"].values, M_corr["lat"].values) and
            np.array_equal(diff_scotese["lon"].values, M_corr["lon"].values)):
        diff_scotese = diff_scotese.interp(
            lat=M_corr["lat"], lon=M_corr["lon"], method="linear",
            kwargs={"fill_value": 0.0})
    M_combined = M_corr + diff_scotese.fillna(0.0)

    # 8. land>0 guard on M_combined (just in case the diff drove land < 0)
    cont_land = (M_orig.values > 0)
    bad = cont_land & (M_combined.values < LAND_MIN_ELEV_M)
    n_bad = int(bad.sum())
    if n_bad > 0:
        if verbose:
            print(f"  ! land>0 guard: clipping {n_bad} cells where dyntopo "
                  f"correction pushed land below sea level "
                  f"(worst {float(M_combined.values[bad].min()):+.0f} m)")
        M_combined.values[bad] = LAND_MIN_ELEV_M

    # 9. write CF-compliant NetCDF
    out_dir.mkdir(parents=True, exist_ok=True)
    out_fp = out_dir / f"{target_age_ma}Ma_corrected_plus_dyntopo_diff_{source}_SW.nc"
    git_hash = git_describe(PALEOTOPO_ROOT)
    source_str = (f"build_dyntopo_diff_correction.py (paleotopo, git {git_hash}); "
                  f"dyntopo from {source} ({_dyntopo_citation(source)})")
    write_cf_grid(
        out_path=out_fp,
        lat=M_corr["lat"].values,
        lon=M_corr["lon"].values,
        fields={
            "M_corrected": (M_corr.values, {
                "units": "m",
                "long_name": "Geochem-corrected paleo-elevation (Zhou et al. 2024 + Scotese & Wright 2018)",
                "standard_name": "surface_altitude",
            }),
            "z_dyntopo_diff": (diff_scotese.values, {
                "units": "m",
                "long_name": f"Per-step dyntopo increment dyntopo({target_age_ma} Ma) − "
                             f"dyntopo({ref_age_ma} Ma), Δt={target_age_ma - ref_age_ma} Myr, "
                             f"plate-frame-computed then cookie-cut + rotated to "
                             f"Scotese paleomag frame at {target_age_ma} Ma",
            }),
            "M_combined": (M_combined.values, {
                "units": "m",
                "long_name": ("M_corrected + dyntopo difference. Cells where "
                              "the dyntopo correction would have pushed "
                              "originally subaerial land (M_orig > 0) below "
                              "sea level are clipped to 0 m — physically "
                              "represents dyntopo-driven coastal inundation "
                              "rather than artificial elevation preservation."),
                "standard_name": "surface_altitude",
            }),
        },
        title=(f"Geochem-corrected Scotese + {source.title()} dyntopo-difference "
               f"paleo-DEM, {target_age_ma} Ma"),
        summary=(f"1° global paleo-DEM with the Zhou et al. (2024) geochem-"
                 f"corrected Scotese & Wright (2018) DEM augmented by the "
                 f"PER-STEP increment of {source.title()} dynamic topography "
                 f"between {ref_age_ma} Ma and {target_age_ma} Ma, evaluated "
                 f"in plate reference frame and cookie-cut + rotated into "
                 f"Scotese paleomag frame at {target_age_ma} Ma."),
        comment=(f"PER-STEP dyntopo coupling: at age t the correction added "
                 f"to M_corrected is Δz(t) = dyntopo(t) − dyntopo(t − Δt), "
                 f"NOT the cumulative-from-present dyntopo(t) − dyntopo(0). "
                 f"At t = 0 the correction is identically zero by "
                 f"construction (M_combined(0) = M_corrected(0) = observed "
                 f"topography); the 0 Ma slice serves as a reference for the "
                 f"t = Δt step only. The subtraction is meaningful only in "
                 f"plate reference frame where each continental block is "
                 f"rigidly attached to itself; the per-step Δz is then "
                 f"cookie-cut by Scotese 2023 continental polygons and "
                 f"rotated to time-t paleomag frame via the Scotese & "
                 f"Wright 2023 rotation model (gplately.Raster.reconstruct). "
                 f"Sea-level guard applied: where the dyntopo correction "
                 f"would have pushed originally subaerial land (M_orig > 0) "
                 f"below sea level, M_combined is clipped to 0 m, "
                 f"representing genuine dyntopo-driven coastal inundation. "
                 f"Cells originally below sea level (M_orig ≤ 0) are "
                 f"unaffected."),
        source=source_str,
        references=("Scotese & Wright 2018 PaleoDEM (Zenodo 5460860); "
                    "Zhou et al. 2024 ESS geochem-assimilated correction; "
                    + _dyntopo_citation(source) + "; "
                    "Scotese & Wright 2023 revised plate model."),
        target_age_ma=target_age_ma,
    )
    if verbose:
        print(f"  wrote {out_fp.name}")
    return out_fp


def _dyntopo_citation(source: str) -> str:
    if source == "young":
        return "Young et al. 2022 G3 dynamic topography (gld428 plate-frame grids)"
    if source == "braz":
        return "Braz et al. 2021 dynamic topography (gmcm9 plate-frame grids)"
    return source


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> int:
    global YOUNG_DIR, BRAZ_DIR
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source", choices=("young", "braz"), default="young",
                    help="Which dyntopo model to use.")
    ap.add_argument("--ages", type=int, nargs="+", default=[5, 15, 35, 55],
                    help="Target ages (Ma) to process.")
    ap.add_argument("--out-dir", type=Path, default=None,
                    help="Override output directory. Default: "
                         "Paleotopo_data_assimilation/data/"
                         "corrected_Scotese_plus_dyntopo_diff_<source>/")
    ap.add_argument("--step-myr", type=int, default=5,
                    help="Δt in Myr — the per-step reference offset. At age "
                         "t the correction added to M_corrected is "
                         "dyntopo(t) − dyntopo(t − Δt), evaluated in plate "
                         "frame. Default 5 Myr.")
    ap.add_argument("--rot-file", type=Path, default=SCOTESE_ROT_FILE)
    ap.add_argument("--polygons", type=Path, default=SCOTESE_POLYGONS)
    ap.add_argument("--young-dir", type=Path, default=YOUNG_DIR,
                    help=f"Override Young 2022 plate-frame grids directory. "
                         f"Default: {YOUNG_DIR}. Can also be set via the "
                         f"PALEOTOPO_DATA_ROOT environment variable.")
    ap.add_argument("--braz-dir", type=Path, default=BRAZ_DIR,
                    help=f"Override Braz 2021 plate-frame grids directory. "
                         f"Default: {BRAZ_DIR}")
    args = ap.parse_args()

    # Mutate module-level globals consulted by the file-lookup helpers.
    YOUNG_DIR = args.young_dir
    BRAZ_DIR  = args.braz_dir

    out_dir = args.out_dir or (PALEOTOPO_ROOT / "data" /
                                f"corrected_Scotese_plus_dyntopo_diff_{args.source}")
    print(f"source         : {args.source}")
    print(f"target ages    : {args.ages}")
    print(f"plate model    : {args.rot_file.name}")
    print(f"polygons       : {args.polygons.name}")
    if args.source == "young":
        print(f"dyntopo dir    : {YOUNG_DIR}")
    else:
        print(f"dyntopo dir    : {BRAZ_DIR}")
    print(f"output dir     : {out_dir}")

    print(f"step Δt        : {args.step_myr} Myr")

    for age in args.ages:
        try:
            process_age(int(age), args.source,
                        rot_file=args.rot_file, polygons_file=args.polygons,
                        out_dir=out_dir, step_myr=args.step_myr)
        except Exception as e:  # noqa: BLE001
            import traceback; traceback.print_exc()
            print(f"FAILED {age} Ma: {type(e).__name__}: {e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
