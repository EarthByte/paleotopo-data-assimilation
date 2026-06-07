"""
=============================================================================
sample_reconstruct_scotese.py  —  Per-slice sample reconstruction (S&W)
=============================================================================

CONCEPTUAL OVERVIEW
    For each target age t, returns the geochem sample paleo-coordinates
    AT THAT EXACT AGE, computed by reconstructing each sample's
    present-day Lat/Lon via the S&W rotation model using a plate ID
    assigned from the Scotese 2008 continental polygons.

    This is the correct way to handle sample positions in the Scotese
    workflow.  The CSV's `PlateID` column was computed against a
    Merdith-family plate model and is NOT consistent with the S&W IDs.
    The S&W topology files only cover the last ~100 Ma so cannot be
    used to look up plate IDs for older slices.

    Instead, we use the Scotese 2008 present-day continental polygons:
    every continental polygon has a `reconstruction_plate_id` that's
    valid for the full age range of the S&W rotation model.  We assign
    each geochem sample the plate ID of whichever continental polygon
    contains its present-day position (point-in-polygon test, done once
    at sample-prep time).  Then for each target age t we reconstruct the
    sample's present-day Lat/Lon using S&W rotations and that assigned
    plate ID.

PUBLIC API
    assign_plate_ids(samples_df, lat_col="Lat", lon_col="Lon",
                     polygons_file=CONTINENTAL_POLYGONS_FILE,
                     rot_model=None)
        Returns an int array of S&W plate IDs, one per sample.  Samples
        in oceanic locations (no continental polygon hit) get 0.

    ScoteseSampleReconstructor()
        Stateful helper.  Call `.assign_plate_ids(df)` once to add an
        `sw_plate_id` column, then `.reconstruct(df, t)` for each
        target age.  Caches rotations across calls.

CACHING
    The FiniteRotation for each (plate_id, target_age) is cached for the
    lifetime of the reconstructor.

PROVENANCE
    Plate IDs are read from `Scotese_Wright_PresentDay_ContinentalPolygons.gpml`
    (or a fallback shapefile of the same name).  The path is set in
    paths_scotese.py.

DEPENDENCIES  pygplates ≥ 1.0, numpy, pandas

USED BY
    assimilate_scotese.py
=============================================================================
"""
from __future__ import annotations
import sys
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import pygplates as pg

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths_scotese import ROT_FILE, CONTINENTAL_POLYGONS_FILE


def load_rotation_model(rot_file: Optional[Path] = None) -> pg.RotationModel:
    return pg.RotationModel(str(rot_file or ROT_FILE))


# ---------------------------------------------------------------------------
# Plate ID assignment via continental polygons
# ---------------------------------------------------------------------------
def assign_plate_ids(samples_df: pd.DataFrame,
                     lat_col: str = "Lat",
                     lon_col: str = "Lon",
                     polygons_file: Optional[Path] = None,
                     rot_model: Optional[pg.RotationModel] = None,
                     ocean_id: int = 0,
                    ) -> np.ndarray:
    """Assign S&W plate IDs to each sample by point-in-polygon against the
    Scotese 2008 present-day continental polygons.

    Parameters
    ----------
    samples_df    : DataFrame containing `lat_col` and `lon_col`
    lat_col       : column name for present-day latitude  (default 'Lat')
    lon_col       : column name for present-day longitude (default 'Lon')
    polygons_file : path to the continental polygons file (defaults to the
                    one configured in paths_scotese.CONTINENTAL_POLYGONS_FILE)
    rot_model     : pre-loaded RotationModel (only required if pyGPlates'
                    PlatePartitioner needs to resolve time-dependent
                    polygons; at t=0 it does not).
    ocean_id      : plate ID assigned to samples not contained in any
                    continental polygon (oceans / shelves / detached).
                    Default 0 means "anchor plate" — those samples will
                    just keep their present-day position when reconstructed.

    Returns
    -------
    plate_ids : np.int32 array, len = len(samples_df)
    """
    pfile = Path(polygons_file or CONTINENTAL_POLYGONS_FILE)
    if not pfile.exists():
        raise FileNotFoundError(
            f"Continental polygons file not found: {pfile}.  "
            "Update CONTINENTAL_POLYGONS_FILE in paths_scotese.py.")
    rot = rot_model or load_rotation_model()

    # PlatePartitioner expects a topological- or static-polygon source plus
    # a rotation model.  For present-day partitioning, time=0.
    partitioner = pg.PlatePartitioner(str(pfile), rot, 0.0)

    out = np.full(len(samples_df), ocean_id, dtype=np.int32)
    for i, row in enumerate(samples_df[[lat_col, lon_col]].itertuples(index=False)):
        la, lo = float(row[0]), float(row[1])
        if not (np.isfinite(la) and np.isfinite(lo)):
            continue
        pt = pg.PointOnSphere((la, lo))
        located = partitioner.partition_point(pt)
        if located is None:
            continue
        try:
            pid = located.get_feature().get_reconstruction_plate_id()
        except Exception:
            continue
        if pid is not None:
            out[i] = int(pid)
    return out


# ---------------------------------------------------------------------------
# Per-slice reconstruction
# ---------------------------------------------------------------------------
class ScoteseSampleReconstructor:
    """Stateful helper bundling plate-ID assignment + per-slice
    reconstruction for the S&W workflow.

    Typical use:
        sr = ScoteseSampleReconstructor()
        df["sw_plate_id"] = sr.assign_plate_ids(df)
        for t in target_ages:
            rlat, rlon = sr.reconstruct(df, t)
    """

    def __init__(self, rot_file: Optional[Path] = None,
                 polygons_file: Optional[Path] = None):
        self.rot_model = load_rotation_model(rot_file)
        self.polygons_file = Path(polygons_file or CONTINENTAL_POLYGONS_FILE)
        self._rot_cache: dict[tuple[int, float], pg.FiniteRotation] = {}

    # ---- plate-ID assignment (one-time, cache on the DataFrame) ----
    def assign_plate_ids(self, samples_df: pd.DataFrame,
                         lat_col: str = "Lat", lon_col: str = "Lon",
                         ocean_id: int = 0) -> np.ndarray:
        return assign_plate_ids(samples_df,
                                lat_col=lat_col, lon_col=lon_col,
                                polygons_file=self.polygons_file,
                                rot_model=self.rot_model,
                                ocean_id=ocean_id)

    # ---- per-slice reconstruction ----
    def reconstruct(self, samples_df: pd.DataFrame, target_age: float,
                    lat_col: str = "Lat", lon_col: str = "Lon",
                    plate_col: str = "sw_plate_id"):
        """Reconstruct each sample to `target_age`.  Samples with
        plate_id == 0 (ocean — not in any continental polygon) are
        returned as NaN; the caller should drop those samples before
        binning into the grid.  Anchor plate is 0 (mantle reference frame).
        """
        n = len(samples_df)
        rlat = np.full(n, np.nan); rlon = np.full(n, np.nan)
        for pid, group in samples_df.groupby(plate_col, dropna=False):
            try:
                pid_int = int(pid)
            except (TypeError, ValueError):
                continue
            if pid_int == 0:
                # leave NaN — these are oceanic samples that should be dropped
                continue
            positions = samples_df.index.get_indexer(group.index)
            key = (pid_int, float(target_age))
            if key not in self._rot_cache:
                try:
                    self._rot_cache[key] = self.rot_model.get_rotation(
                        float(target_age), pid_int, anchor_plate_id=0)
                except Exception:
                    self._rot_cache[key] = None
            rot = self._rot_cache[key]
            if rot is None:
                continue
            for pos, la, lo in zip(positions,
                                   group[lat_col].to_numpy(),
                                   group[lon_col].to_numpy()):
                try:
                    p_t = rot * pg.PointOnSphere((float(la), float(lo)))
                    la_t, lo_t = p_t.to_lat_lon()
                    rlat[pos] = la_t; rlon[pos] = lo_t
                except Exception:
                    pass
        return rlat, rlon


if __name__ == "__main__":
    sr = ScoteseSampleReconstructor()
    df = pd.DataFrame({
        "Lat":[40.0, -23.0, 30.0],     # N. America Plains, Brazil interior, India
        "Lon":[-100.0, -50.0, 78.0],
    })
    pids = sr.assign_plate_ids(df)
    df["sw_plate_id"] = pids
    print("assigned plate IDs:", pids)
    for t in [0, 50, 150, 300, 500]:
        rla, rlo = sr.reconstruct(df, t)
        print(f"t={t:4d} Ma →",
              [(round(rla[i],1), round(rlo[i],1)) for i in range(len(df))])
