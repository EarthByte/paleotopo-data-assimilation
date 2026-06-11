#!/usr/bin/env python3
"""
Rotate Young et al. (2022) dynamic-topography grids from Muller2022 NNR
to the Scotese paleomag frame at each grid's own age (gplately).

Two source modes are supported:

  --src-5myr <dir>    *Production path.* Native 5-Myr release
                      (gld428DT_MFgrid<age>_corrected.nc, 0..250 Ma).
                      Each grid is rotated at its own age and written
                      directly to the destination with no time
                      interpolation. Sanity-check videos confirmed that
                      this is the right pipeline; 5-Myr-native rotated
                      polygons line up with reconstructed Scotese
                      continental polygons at every age.

  --src-20myr <dir>   *Deprecated.* Native 20-Myr release
                      (gld428DT_MantleFrameGrid<age>.grd, 0..580+ Ma).
                      Rotated grids are then linearly time-interpolated
                      per-cell to a 5-Myr cadence (Stage B). Initially
                      attractive because the 20-Myr release is smoother,
                      but linear blending between Scotese-frame rotated
                      grids at consecutive 20-Myr nodes is *not* a valid
                      paleo-frame state at the intermediate ages:
                      Scotese-frame plate motion within each 20-Myr
                      window has nowhere to go in the blend, so the
                      blended grid is geographically off at e.g. 5, 15,
                      25 Ma. Kept here only for reference / sanity
                      checks. Do not use as the production path.

Frame change rationale (same as porphyry_grids_scotese.py):
  source : Muller2022/optimisation/no_net_rotation_model.rot @ anchor 0
  target : Scotese_Wright_PlateModel_2023.rot @ anchor 0 (paleomag native)

Valid for the last 300 Ma. The 0–250 Ma window is the Scotese ↔ Merdith
plate-motion convergence window where the 5-Myr-native release is used.
The 255–300 Ma window is the Pangean stability interval where plates are
not moving much, so the 20-Myr-release + per-cell linear time
interpolation between bracketing 20-Myr nodes is a defensible
approximation (no significant Scotese-frame plate motion within each
20-Myr window). Beyond 300 Ma the dyntopo variant is not produced.

Output
------

`<dst>/dyntopo_scotese_<age>Ma.nc` — one NetCDF per 5-Myr age slice in
[--age-min, --age-max] (capped to 250). Variable `z` is dynamic
topography in metres. Grid: 0.5° resolution (361 × 721),
-90°…+90° lat × -180°…+180° lon.

Usage
-----

  # production path: 5-Myr native source, rotate at each age
  python scripts_Scotese/rotate_young_dyntopo_to_scotese.py \\
      --src-5myr data/Young2022_gld428_grids_5Myr \\
      --dst data/Young2022_gld428_grids_5Myr_scotese_frame \\
      --muller2022-dir data/Muller_etal_2022_SE_1Ga_Opt_PlateMotionModel_v1.2.4 \\
      --scotese-dir data/Scotese_Wright_revised_plate_model \\
      --age-min 0 --age-max 250 --cadence-myr 5
"""

from __future__ import annotations

import argparse
import re
import sys
import tempfile
import time
from pathlib import Path

import numpy as np


SRC_20MYR_FNAME_RE = re.compile(r"^gld428DT_MantleFrameGrid(\d+)\.grd$")
SRC_20MYR_FNAME_FMT = "gld428DT_MantleFrameGrid{age}.grd"
SRC_5MYR_FNAME_RE = re.compile(r"^gld428DT_MFgrid(\d+)_corrected\.nc$")
SRC_5MYR_FNAME_FMT = "gld428DT_MFgrid{age}_corrected.nc"

DST_FNAME_FMT = "dyntopo_scotese_{age}Ma.nc"
ROTATED_INTERMEDIATE_FMT = "dyntopo_scotese_native_{age}Ma.nc"

MULLER2022_NNR_REL = "optimisation/no_net_rotation_model.rot"
SCOTESE_ROT_FILE = "Scotese_Wright_PlateModel_2023.rot"


# ----------------------------------------------------------------------------
# Stage A: rotate one 20-Myr grid (worker-safe for joblib)
# ----------------------------------------------------------------------------

def _rotate_one_grid(input_path: str, output_path: str,
                     reconstruction_time: float,
                     from_rotation_files: list[str],
                     to_rotation_files: list[str],
                     from_anchor: int, to_anchor: int,
                     grid_spacing_degrees: float) -> str:
    import contextlib
    import io
    import warnings

    import gplately
    import gplately.pygplates as pygplates

    with contextlib.redirect_stdout(io.StringIO()), \
         contextlib.redirect_stderr(io.StringIO()), \
         warnings.catch_warnings():
        warnings.simplefilter("ignore")
        src_model = pygplates.RotationModel([str(p) for p in from_rotation_files])
        tgt_model = pygplates.RotationModel([str(p) for p in to_rotation_files])
        raster = gplately.grids.Raster(input_path, time=float(reconstruction_time))
        raster.rotate_reference_frames(
            grid_spacing_degrees=grid_spacing_degrees,
            reconstruction_time=float(reconstruction_time),
            from_rotation_features_or_model=src_model,
            to_rotation_features_or_model=tgt_model,
            from_rotation_reference_plate=int(from_anchor),
            to_rotation_reference_plate=int(to_anchor),
            output_name=output_path,
        )
    return output_path


# ----------------------------------------------------------------------------
# Stage B: time-interpolation in Scotese frame
# ----------------------------------------------------------------------------

def _load_rotated_grid(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read a rotated NetCDF and return (lat, lon, z) arrays."""
    from netCDF4 import Dataset
    with Dataset(path, "r") as ds:
        lat = np.asarray(ds["lat"][:], dtype=np.float64)
        lon = np.asarray(ds["lon"][:], dtype=np.float64)
        z   = np.asarray(ds["z"][:],   dtype=np.float64)
    return lat, lon, z


def _save_grid(path: Path, lat: np.ndarray, lon: np.ndarray, z: np.ndarray,
               age_ma: int, source_ages: tuple[int, int],
               source_release: str) -> None:
    """Write a NetCDF with z(lat, lon) in matching schema to
    rotated source grids, plus metadata documenting source + method.

    source_release : "20myr" or "5myr_native". Controls the method
                     attribute and the name of the source-ages attribute.
    """
    from netCDF4 import Dataset
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with Dataset(tmp, "w") as ds:
        ds.createDimension("lat", lat.size)
        ds.createDimension("lon", lon.size)
        vlat = ds.createVariable("lat", "f8", ("lat",))
        vlon = ds.createVariable("lon", "f8", ("lon",))
        vz   = ds.createVariable("z",   "f4", ("lat", "lon"))
        vlat[:] = lat
        vlon[:] = lon
        vz[:, :] = z.astype(np.float32)
        ds.age_ma = int(age_ma)
        if source_release == "20myr":
            ds.source_20myr_ages = list(source_ages)
            ds.method = (
                "Young 2022 dyntopo (20-Myr release) @ Muller2022 NNR → "
                "Scotese paleomag via gplately Raster.rotate_reference_frames; "
                "linearly interpolated between the two bracketing 20-Myr "
                "rotated grids."
            )
        elif source_release == "5myr_native":
            ds.source_5myr_age = int(source_ages[0])
            ds.method = (
                "Young 2022 dyntopo (5-Myr native release, *_corrected.nc) @ "
                "Muller2022 NNR → Scotese paleomag via gplately "
                "Raster.rotate_reference_frames. No time interpolation."
            )
        else:
            raise ValueError(f"unknown source_release {source_release!r}")
    tmp.replace(path)


def _interpolate(age_target: int, native_ages: list[int],
                 rotated_dir: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, tuple[int, int]]:
    """Return (lat, lon, z, (a_lo, a_hi)) for a target age, by linear
    per-cell interpolation between bracketing rotated grids."""
    if age_target in native_ages:
        lat, lon, z = _load_rotated_grid(
            rotated_dir / ROTATED_INTERMEDIATE_FMT.format(age=age_target),
        )
        return lat, lon, z, (age_target, age_target)

    # Find bracketing 20-Myr ages
    arr = np.array(native_ages)
    below = arr[arr < age_target]
    above = arr[arr > age_target]
    if below.size == 0 or above.size == 0:
        raise ValueError(
            f"Cannot bracket {age_target} Ma; native ages range "
            f"{arr.min()}..{arr.max()}"
        )
    a_lo = int(below.max())
    a_hi = int(above.min())

    lat_lo, lon_lo, z_lo = _load_rotated_grid(
        rotated_dir / ROTATED_INTERMEDIATE_FMT.format(age=a_lo),
    )
    _, _, z_hi = _load_rotated_grid(
        rotated_dir / ROTATED_INTERMEDIATE_FMT.format(age=a_hi),
    )
    w = (age_target - a_lo) / (a_hi - a_lo)
    z_interp = (1.0 - w) * z_lo + w * z_hi
    return lat_lo, lon_lo, z_interp, (a_lo, a_hi)


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

def _run_20myr_path(args, nnr_rot: Path, sco_rot: Path) -> int:
    """Production path: rotate 20-Myr release, then interpolate to cadence."""
    # Pick source 20-Myr ages that bracket the target age range.
    # Need one below age_min and one above age_max to interpolate the edges.
    src_ages_avail: list[int] = sorted(
        int(SRC_20MYR_FNAME_RE.match(p.name).group(1))
        for p in args.src_20myr.iterdir()
        if SRC_20MYR_FNAME_RE.match(p.name)
    )
    if not src_ages_avail:
        raise SystemExit(f"No gld428DT_MantleFrameGrid<age>.grd under {args.src_20myr}")
    print(f"20-Myr source ages available ({len(src_ages_avail)}): "
          f"{src_ages_avail[0]}–{src_ages_avail[-1]} Ma")

    # Native 20-Myr ages we actually need to rotate: those between age_min and
    # age_max, plus one beyond each end for edge-interpolation.
    need_native = [
        a for a in src_ages_avail
        if (args.age_min - 20) <= a <= (min(args.age_max, 540) + 20)
    ]
    print(f"Native 20-Myr ages to rotate ({len(need_native)}): "
          f"{need_native[0]}–{need_native[-1]} Ma")

    # ----- Stage A: rotate native 20-Myr grids -----
    # Use a temp dir for intermediates (rotated 20-Myr grids in Scotese frame)
    if args.keep_intermediate:
        intermediate_dir = args.dst / "rotated_20myr"
        intermediate_dir.mkdir(parents=True, exist_ok=True)
        cleanup = False
    else:
        intermediate_dir = Path(tempfile.mkdtemp(prefix="rotated_20myr_"))
        cleanup = True

    todo_rotate: list[int] = []
    for age in need_native:
        out = intermediate_dir / ROTATED_INTERMEDIATE_FMT.format(age=age)
        if not out.exists():
            todo_rotate.append(age)
    if todo_rotate:
        print(f"\n[Stage A] rotating {len(todo_rotate)} native 20-Myr grids "
              f"({len(need_native) - len(todo_rotate)} already done)")
        # Prime models in main process
        import gplately.pygplates as pygplates
        _ = pygplates.RotationModel(str(nnr_rot))
        _ = pygplates.RotationModel(str(sco_rot))

        from joblib import Parallel, delayed
        t0 = time.time()
        Parallel(n_jobs=args.n_jobs, backend="loky", verbose=1)(
            delayed(_rotate_one_grid)(
                input_path=str(args.src_20myr / SRC_20MYR_FNAME_FMT.format(age=age)),
                output_path=str(intermediate_dir / ROTATED_INTERMEDIATE_FMT.format(age=age)),
                reconstruction_time=float(age),
                from_rotation_files=[str(nnr_rot)],
                to_rotation_files=[str(sco_rot)],
                from_anchor=0,
                to_anchor=0,
                grid_spacing_degrees=args.grid_spacing_degrees,
            )
            for age in todo_rotate
        )
        dt = time.time() - t0
        n_native_ok = sum(1 for a in todo_rotate
                          if (intermediate_dir /
                              ROTATED_INTERMEDIATE_FMT.format(age=a)).exists())
        print(f"  [Stage A] {n_native_ok}/{len(todo_rotate)} rotated in {dt/60:.1f} min")
    else:
        print("\n[Stage A] all native 20-Myr rotations already done")

    # ----- Stage B: time-interpolate to cadence -----
    target_ages = list(range(args.age_min,
                             min(args.age_max, 540) + 1,
                             args.cadence_myr))
    print(f"\n[Stage B] interpolating to {len(target_ages)} target ages "
          f"({target_ages[0]}–{target_ages[-1]} Ma at {args.cadence_myr} Myr cadence)")

    n_interp = n_skip = 0
    for age_target in target_ages:
        out_path = args.dst / DST_FNAME_FMT.format(age=age_target)
        if out_path.exists():
            n_skip += 1
            continue
        try:
            lat, lon, z, src_pair = _interpolate(age_target, need_native, intermediate_dir)
        except Exception as e:  # noqa: BLE001
            print(f"  age {age_target} Ma: FAIL ({e})")
            continue
        _save_grid(out_path, lat, lon, z,
                   age_ma=age_target, source_ages=src_pair,
                   source_release="20myr")
        n_interp += 1
        if age_target % 50 == 0:
            print(f"  age {age_target} Ma: interpolated from {src_pair[0]}-{src_pair[1]} Ma  "
                  f"z range [{z.min():.0f}, {z.max():.0f}] m")
    print(f"\n[Stage B] {n_interp} interpolated, {n_skip} skipped")

    # ----- Cleanup intermediates -----
    if cleanup:
        import shutil
        shutil.rmtree(intermediate_dir, ignore_errors=True)
        print(f"\n[cleanup] removed temp intermediates at {intermediate_dir}")
    else:
        print(f"\n[cleanup] kept intermediates at {intermediate_dir} (--keep-intermediate)")

    print(f"\nDone. Final {args.cadence_myr}-Myr-cadence grids (from 20-Myr "
          f"source + interpolation): {args.dst}")
    return 0


def _run_5myr_native_path(args, nnr_rot: Path, sco_rot: Path) -> int:
    """Sanity-check path: rotate native 5-Myr release, no interpolation."""
    src_ages_avail: list[int] = sorted(
        int(SRC_5MYR_FNAME_RE.match(p.name).group(1))
        for p in args.src_5myr.iterdir()
        if SRC_5MYR_FNAME_RE.match(p.name)
    )
    if not src_ages_avail:
        raise SystemExit(f"No gld428DT_MFgrid<age>_corrected.nc under {args.src_5myr}")
    print(f"5-Myr native source ages available ({len(src_ages_avail)}): "
          f"{src_ages_avail[0]}–{src_ages_avail[-1]} Ma")

    target_ages = [
        a for a in range(args.age_min,
                         min(args.age_max, 540) + 1,
                         args.cadence_myr)
        if a in src_ages_avail
    ]
    missing = [a for a in range(args.age_min,
                                min(args.age_max, 540) + 1,
                                args.cadence_myr)
               if a not in src_ages_avail]
    if missing:
        print(f"WARNING: requested ages with no 5-Myr native source: {missing}")
    print(f"Native 5-Myr ages to rotate ({len(target_ages)}): "
          f"{target_ages[0]}–{target_ages[-1]} Ma at {args.cadence_myr} Myr cadence")

    # Rotate each native 5-Myr grid directly into a temp intermediate, then
    # re-save through _save_grid so the destination has our standard schema
    # and metadata.
    if args.keep_intermediate:
        intermediate_dir = args.dst / "rotated_5myr_native"
        intermediate_dir.mkdir(parents=True, exist_ok=True)
        cleanup = False
    else:
        intermediate_dir = Path(tempfile.mkdtemp(prefix="rotated_5myr_native_"))
        cleanup = True

    todo_rotate: list[int] = []
    for age in target_ages:
        out = args.dst / DST_FNAME_FMT.format(age=age)
        if not out.exists():
            todo_rotate.append(age)
    if not todo_rotate:
        print("\nAll target ages already present in dst; nothing to do.")
        if cleanup:
            import shutil
            shutil.rmtree(intermediate_dir, ignore_errors=True)
        return 0

    print(f"\n[rotate] {len(todo_rotate)} native 5-Myr grids "
          f"({len(target_ages) - len(todo_rotate)} already done)")
    # Prime models in main process
    import gplately.pygplates as pygplates
    _ = pygplates.RotationModel(str(nnr_rot))
    _ = pygplates.RotationModel(str(sco_rot))

    from joblib import Parallel, delayed
    t0 = time.time()
    Parallel(n_jobs=args.n_jobs, backend="loky", verbose=1)(
        delayed(_rotate_one_grid)(
            input_path=str(args.src_5myr / SRC_5MYR_FNAME_FMT.format(age=age)),
            output_path=str(intermediate_dir / ROTATED_INTERMEDIATE_FMT.format(age=age)),
            reconstruction_time=float(age),
            from_rotation_files=[str(nnr_rot)],
            to_rotation_files=[str(sco_rot)],
            from_anchor=0,
            to_anchor=0,
            grid_spacing_degrees=args.grid_spacing_degrees,
        )
        for age in todo_rotate
    )
    dt = time.time() - t0
    n_ok = sum(1 for a in todo_rotate
               if (intermediate_dir / ROTATED_INTERMEDIATE_FMT.format(age=a)).exists())
    print(f"  [rotate] {n_ok}/{len(todo_rotate)} rotated in {dt/60:.1f} min")

    # ----- Re-save with our schema/metadata -----
    n_saved = 0
    for age in todo_rotate:
        src = intermediate_dir / ROTATED_INTERMEDIATE_FMT.format(age=age)
        if not src.exists():
            print(f"  age {age} Ma: missing rotated intermediate, skipping")
            continue
        lat, lon, z = _load_rotated_grid(src)
        out_path = args.dst / DST_FNAME_FMT.format(age=age)
        _save_grid(out_path, lat, lon, z,
                   age_ma=age, source_ages=(age, age),
                   source_release="5myr_native")
        n_saved += 1
        if age % 50 == 0:
            print(f"  age {age} Ma: z range [{z.min():.0f}, {z.max():.0f}] m")
    print(f"\n[save] {n_saved} written to {args.dst}")

    # ----- Cleanup intermediates -----
    if cleanup:
        import shutil
        shutil.rmtree(intermediate_dir, ignore_errors=True)
        print(f"\n[cleanup] removed temp intermediates at {intermediate_dir}")
    else:
        print(f"\n[cleanup] kept intermediates at {intermediate_dir} (--keep-intermediate)")

    print(f"\nDone. 5-Myr native rotated grids: {args.dst}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Rotate Young 2022 dyntopo grids from Muller2022 NNR to "
                    "Scotese paleomag frame. With --src-20myr the rotated "
                    "grids are time-interpolated to --cadence-myr (production "
                    "path). With --src-5myr the native 5-Myr release is "
                    "rotated directly with no interpolation (sanity-check)."
    )
    src_grp = ap.add_mutually_exclusive_group(required=True)
    src_grp.add_argument("--src-20myr", type=Path,
                         help="Directory of gld428DT_MantleFrameGrid<age>.grd files")
    src_grp.add_argument("--src-5myr", type=Path,
                         help="Directory of gld428DT_MFgrid<age>_corrected.nc files")
    ap.add_argument("--dst", type=Path, required=True,
                    help="Output directory for rotated NetCDFs")
    ap.add_argument("--muller2022-dir", type=Path, required=True)
    ap.add_argument("--scotese-dir", type=Path, required=True)
    ap.add_argument("--age-min", type=int, default=0)
    ap.add_argument("--age-max", type=int, default=250,
                    help="Cap (Scotese↔Merdith plate-motion convergence window)")
    ap.add_argument("--cadence-myr", type=int, default=5,
                    help="Output cadence in Ma; default 5")
    ap.add_argument("--grid-spacing-degrees", type=float, default=0.5,
                    help="Rotation output grid resolution; default 0.5° = native")
    ap.add_argument("--n-jobs", type=int, default=-2,
                    help="joblib n_jobs for rotations (default -2)")
    ap.add_argument("--keep-intermediate", action="store_true",
                    help="Don't delete rotated intermediates after the run")
    args = ap.parse_args()

    src_dir = args.src_20myr if args.src_20myr else args.src_5myr
    if not src_dir.is_dir():
        raise SystemExit(f"source directory not found: {src_dir}")
    args.dst.mkdir(parents=True, exist_ok=True)

    nnr_rot = args.muller2022_dir / MULLER2022_NNR_REL
    sco_rot = args.scotese_dir / SCOTESE_ROT_FILE
    for p in (nnr_rot, sco_rot):
        if not p.exists():
            raise SystemExit(f"Missing rotation file: {p}")

    if args.src_20myr:
        return _run_20myr_path(args, nnr_rot, sco_rot)
    else:
        return _run_5myr_native_path(args, nnr_rot, sco_rot)


if __name__ == "__main__":
    sys.exit(main())
