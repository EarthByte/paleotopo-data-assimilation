"""
=============================================================================
plate_model_utils_scotese.py  —  pyGPlates helpers (S&W workflow)
=============================================================================

Same purpose as plate_model_utils.py in the Merdith workflow:  turn the
plate model into 1°-grid rasters at any geological age t.

Three public functions:

  cob_mask(t, lat, lon, M)     -> bool[n_lat, n_lon]
        Continental footprint at age t.  S&W does NOT ship a separate
        continent-ocean-boundary shapefile, so the mask is derived from
        the PaleoDEM itself: cells with elevation > COB_PROXY_M are
        treated as continental.  This is robust and avoids the
        terrane-classification ambiguity of using the plate polygons.

  subduction_zones(t)          -> list[(lat,lon) point arrays]
        Active subduction-zone segments at age t resolved from the S&W
        topology files.

  province_grid(t, lat, lon, cob_mask, declustered_samples, M)
        int8[n_lat, n_lon] tectonic-province raster.  Logic is identical
        to the Merdith version (margin → arc → geomorphic surrogate →
        sample-data overlay), but uses S&W-derived continent mask and
        S&W-derived subduction zones.

CONFIGURATION  (top-of-file constants)
    ARC_DIST_KM     = 450 km  → "Continental Arc" province threshold
    MARGIN_DIST_KM  = 200 km  → "Continental Margin" province threshold
    COB_PROXY_M     = -1500 m → continental footprint threshold for the
                                DEM-based COB mask (Scotese & Wright)
    Geomorphic surrogate thresholds:
        ELEV_OROGEN_M    = 1200    cell elevation > this → Orogen
        ROUGH_OROGEN_M   = 250     local σ(elev) > this → Orogen
        ELEV_SHIELD_HI_M = 600     cell elevation < this AND
        ROUGH_SHIELD_LO_M= 100     local σ < this → Shield

DEPENDENCIES
    pygplates 1.0+, numpy, scipy, matplotlib

USED BY
    assimilate_scotese.py, render_v2_scotese.py, render_videos_pygmt_scotese.py
=============================================================================
"""
from __future__ import annotations
import os, math, sys
from pathlib import Path
import numpy as np
import pygplates as pg
from matplotlib.path import Path as MplPath

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths_scotese import (
    PLATE_DIR, ROT_FILE, BOUNDARIES_FILE, POLYGONS_FILE,
    CONTINENTAL_POLYGONS_FILE,
)

TOPOLOGY_FILES = [BOUNDARIES_FILE, POLYGONS_FILE]

ARC_DIST_KM = 450.0
MARGIN_DIST_KM = 200.0

# Fallback DEM-based COB threshold (only used if the continental polygons
# file is missing).  −1500 m captures continental shelves robustly.
COB_PROXY_M = -1500.0

PROV_LIST = ["Continental Arc", "Orogen", "Continental Margin",
             "Island Arc", "Extended Crust", "Basin", "Platform", "Shield", "Other"]
PROV_INDEX = {p: i for i, p in enumerate(PROV_LIST)}


# ---------------------------------------------------------------------------
# Rotation model (lazy load + cache)
# ---------------------------------------------------------------------------
_rot_cache = None

def get_rotation_model():
    global _rot_cache
    if _rot_cache is None:
        _rot_cache = pg.RotationModel(str(ROT_FILE))
    return _rot_cache


# ---------------------------------------------------------------------------
# Continent mask
#   Primary path: rasterise the Scotese 2008 continental polygons
#                 reconstructed via the S&W rotation model to age t.
#   Fallback:     DEM threshold (z > COB_PROXY_M, default −1500 m) — used
#                 only if the continental polygons file is missing.
# ---------------------------------------------------------------------------
def _polygon_to_latlon(poly_on_sphere):
    pts = [p.to_lat_lon() for p in poly_on_sphere]
    return np.array(pts)  # (N, 2): (lat, lon)


def _split_dateline(latlon):
    """Split a polygon at any dateline crossing so matplotlib Path works.
    Returns list of (lat, lon) sub-polygons."""
    lat_arr = latlon[:, 0]; lon_arr = latlon[:, 1].copy()
    if len(lon_arr) < 3:
        return []
    dlon = np.diff(lon_arr)
    if np.max(np.abs(dlon)) < 180:
        return [latlon]
    lon360 = np.where(lon_arr < 0, lon_arr + 360, lon_arr)
    if np.max(np.abs(np.diff(lon360))) < 180:
        return [np.column_stack([lat_arr, lon360 - 180])]
    idx = np.argmax(np.abs(dlon)) + 1
    a = latlon[:idx]; b = latlon[idx:]
    return [a, b] if len(a) >= 3 and len(b) >= 3 else [latlon]


def _rasterize_polygons(reconstructed, lat, lon):
    """True where any reconstructed PolygonOnSphere covers a grid cell centre."""
    nlat, nlon = len(lat), len(lon)
    LAT2D, LON2D = np.meshgrid(lat, lon, indexing="ij")
    pts = np.column_stack([LON2D.ravel(), LAT2D.ravel()])
    mask = np.zeros((nlat, nlon), dtype=bool)
    for rf in reconstructed:
        geom = rf.get_reconstructed_geometry()
        if not isinstance(geom, pg.PolygonOnSphere):
            continue
        latlon = _polygon_to_latlon(geom)
        for sub in _split_dateline(latlon):
            if len(sub) < 3: continue
            path = MplPath(sub[:, [1, 0]])
            inside = path.contains_points(pts).reshape(nlat, nlon)
            mask |= inside
    return mask


def cob_mask(t: float, lat: np.ndarray, lon: np.ndarray,
             M: np.ndarray = None, threshold_m: float = COB_PROXY_M,
             use_polygons: bool = True) -> np.ndarray:
    """Continental footprint at age t.

    By default the Scotese 2008 continental polygons are reconstructed via
    the S&W rotation model and rasterised onto the lat/lon grid.  If the
    polygons file is missing (or `use_polygons=False`), falls back to a
    DEM threshold (`M > threshold_m`).

    Parameters
    ----------
    t           : age in Ma
    lat, lon    : 1-D grid axes
    M           : (optional) elevation grid, only consulted in the fallback
    threshold_m : DEM-fallback elevation cutoff (default −1500 m)
    use_polygons: if True, try polygons first; if False go straight to DEM
    """
    if use_polygons and CONTINENTAL_POLYGONS_FILE.exists():
        rot = get_rotation_model()
        out = []
        pg.reconstruct(str(CONTINENTAL_POLYGONS_FILE), rot, out, float(t))
        return _rasterize_polygons(out, lat, lon)
    if M is None:
        raise ValueError("DEM elevation grid M required for fallback mask")
    return M > threshold_m


# ---------------------------------------------------------------------------
# Subduction zones
# ---------------------------------------------------------------------------
def _resolved_topologies(t):
    rot = get_rotation_model()
    resolved = []; shared = []
    pg.resolve_topologies([str(f) for f in TOPOLOGY_FILES], rot, resolved, t, shared)
    return resolved, shared


def subduction_zones(t):
    """Return list of (lat,lon) arrays for each subduction-zone segment at t."""
    _, shared = _resolved_topologies(t)
    out = []
    for sb in shared:
        feat = sb.get_feature()
        if feat.get_feature_type() != pg.FeatureType.gpml_subduction_zone:
            continue
        for sub in sb.get_shared_sub_segments():
            geom = sub.get_resolved_geometry()
            if geom is None:
                continue
            try:
                pts = np.array([p.to_lat_lon() for p in geom])
            except TypeError:
                continue
            if len(pts) >= 2:
                out.append(pts)
    return out


def plate_polygons(t):
    """Closed polygon outlines of every RESOLVED plate-topology polygon
    at time t.

    Returns
    -------
    list of (N, 2) float arrays, one per resolved topology
        Each ring is in [lat, lon] order (degrees), last vertex equal to
        first so an ``ax.plot`` draws the full closed loop.

    Mirrors ``plate_model_utils.plate_polygons`` in the Merdith
    workflow.  Used by the renderer overlays as a guaranteed-closed
    backstop for the SZ-line + teeth overlay — filtering by
    feature type alone leaves visible gaps in the plate-boundary
    rendering where unclassified/unlabelled topology features are
    present.  In the S&W topology set this matters less than for
    Merdith (S&W resolves topologies only for ages <= ~100 Ma) but the
    same fix applies.

    Returns an empty list at ages where the S&W topology files do not
    resolve (older than ~100 Ma); callers should handle the empty case
    by skipping the overlay rather than raising.
    """
    resolved, _ = _resolved_topologies(t)
    rings = []
    for rb in resolved:
        geom = rb.get_resolved_geometry()
        if not isinstance(geom, pg.PolygonOnSphere):
            continue
        try:
            latlon = np.array([p.to_lat_lon() for p in geom])
        except TypeError:
            continue
        if len(latlon) < 3:
            continue
        if not np.allclose(latlon[0], latlon[-1]):
            latlon = np.vstack([latlon, latlon[0:1]])
        rings.append(latlon)
    return rings


def topology_lines(t):
    """All plate-boundary topology line segments at time t, categorised
    for the renderer overlay.  Mirrors ``plate_model_utils.topology_lines``
    in the Merdith workflow.

    Returns a dict
      "polygons":    list of closed-ring (N, 2) [lat, lon] arrays — one
                     per resolved plate topology.  Drawn as a thin
                     backstop by the renderers so boundaries close
                     regardless of per-segment feature-type labelling.
      "subduction":  list of (pts, polarity) tuples — for SZ teeth.
      "other":       list of pts arrays for every non-subduction
                     boundary segment (mid-ocean ridges, transforms,
                     and unclassified boundaries lumped here).

    At ages beyond the S&W topology coverage (~100 Ma), all three
    lists come back empty.
    """
    resolved, shared = _resolved_topologies(t)
    out = {"polygons": [], "subduction": [], "other": []}
    pol_prop = pg.PropertyName.create_gpml("subductionPolarity")

    for rb in resolved:
        geom = rb.get_resolved_geometry()
        if not isinstance(geom, pg.PolygonOnSphere):
            continue
        try:
            latlon = np.array([p.to_lat_lon() for p in geom])
        except TypeError:
            continue
        if len(latlon) < 3:
            continue
        if not np.allclose(latlon[0], latlon[-1]):
            latlon = np.vstack([latlon, latlon[0:1]])
        out["polygons"].append(latlon)

    for sb in shared:
        feat_type = sb.get_feature().get_feature_type()
        is_sub = (feat_type == pg.FeatureType.gpml_subduction_zone)
        for sub in sb.get_shared_sub_segments():
            geom = sub.get_resolved_geometry()
            if geom is None:
                continue
            try:
                pts = np.array([p.to_lat_lon() for p in geom])
            except TypeError:
                continue
            if len(pts) < 2:
                continue
            if is_sub:
                try:
                    pol = sub.get_feature().get_enumeration(pol_prop, "Unknown")
                except Exception:
                    pol = "Unknown"
                out["subduction"].append((pts, pol))
            else:
                out["other"].append(pts)
    return out


def distance_to_subduction(t, lat, lon, sample_step=4):
    """Min distance (km) from each grid cell to the nearest active SZ
    segment at time t.  Subsamples grid + SZ vertices for speed; the
    Gaussian arc kernel is much wider than the residual error."""
    sz_lines = subduction_zones(t)
    nlat, nlon = len(lat), len(lon)
    if not sz_lines:
        return np.full((nlat, nlon), np.inf)
    all_pts = np.concatenate(sz_lines, axis=0)
    step_pt = max(1, len(all_pts) // 400)
    all_pts = all_pts[::step_pt]

    R = 6371.0
    LAT2D, LON2D = np.meshgrid(lat[::sample_step], lon[::sample_step], indexing="ij")
    lat_g = np.radians(LAT2D)[..., None]
    lon_g = np.radians(LON2D)[..., None]
    lat_p = np.radians(all_pts[:, 0])[None, None, :]
    lon_p = np.radians(all_pts[:, 1])[None, None, :]
    dlat = lat_p - lat_g
    dlon = lon_p - lon_g
    a = np.sin(dlat/2)**2 + np.cos(lat_g)*np.cos(lat_p)*np.sin(dlon/2)**2
    d_sub = 2 * R * np.arcsin(np.sqrt(np.clip(a, 0, 1)))
    d_min_sub = d_sub.min(axis=-1)
    # upsample back to full grid via nearest-neighbour replication
    d_full = np.repeat(np.repeat(d_min_sub, sample_step, axis=0), sample_step, axis=1)
    return d_full[:nlat, :nlon]


# ---------------------------------------------------------------------------
# Province grid (geomorphic surrogate + plate-model arc + sample overlay)
# ---------------------------------------------------------------------------
def local_std(M, mask, k=2):
    """Fast local stdev over (2k+1)x(2k+1) window, NaN where mask=False."""
    from scipy.ndimage import uniform_filter
    M_filled = np.where(mask, M.astype(float), 0.0)
    cnt = uniform_filter(mask.astype(float), size=2*k+1, mode="nearest")
    mean = uniform_filter(M_filled, size=2*k+1, mode="nearest") / np.maximum(cnt, 1e-9)
    sq = uniform_filter(M_filled**2, size=2*k+1, mode="nearest") / np.maximum(cnt, 1e-9)
    var = sq - mean**2
    out = np.sqrt(np.maximum(var, 0.0))
    out[~mask] = np.nan
    return out


def province_grid(t, lat, lon, cont_mask, dec, M,
                  arc_dist_km=ARC_DIST_KM, margin_dist_km=MARGIN_DIST_KM,
                  elev_orogen=1200.0, rough_orogen=250.0,
                  elev_shield_hi=600.0, rough_shield_lo=100.0):
    """S&W plate-model + geomorphic-fallback province classifier."""
    from scipy.ndimage import distance_transform_edt
    nlat, nlon = M.shape
    P_grid = np.full((nlat, nlon), PROV_INDEX["Other"], dtype=np.int8)
    P_grid[~cont_mask] = -1

    # Margin: continental cells within ~margin_dist_km of the COB outer edge
    if cont_mask.any():
        d_to_ocean_cells = distance_transform_edt(cont_mask)
        d_to_ocean_km = d_to_ocean_cells * 111.32
        margin_mask = cont_mask & (d_to_ocean_km < margin_dist_km)
        P_grid[margin_mask] = PROV_INDEX["Continental Margin"]

    # Arc: continental cells within arc_dist_km of an active SZ
    d_sz = distance_to_subduction(t, lat, lon, sample_step=4)
    arc_mask = cont_mask & (d_sz < arc_dist_km)
    P_grid[arc_mask] = PROV_INDEX["Continental Arc"]

    # Geomorphic surrogate inside the unlabelled continental interior
    rough = local_std(M, cont_mask, k=2)
    interior = cont_mask & (P_grid == PROV_INDEX["Other"])
    orogen_mask = interior & ((M > elev_orogen) | ((rough > rough_orogen) & (M > 200)))
    P_grid[orogen_mask] = PROV_INDEX["Orogen"]
    shield_mask = interior & (M < elev_shield_hi) & (M > 0) & (rough < rough_shield_lo)
    P_grid[shield_mask] = PROV_INDEX["Shield"]
    platform_mask = interior & (P_grid == PROV_INDEX["Other"]) & (M > 0)
    P_grid[platform_mask] = PROV_INDEX["Platform"]

    # Sample-province overlay in small halos around their paleo locations.
    paint_priority = ["Shield", "Platform", "Extended Crust", "Basin",
                      "Island Arc", "Continental Arc", "Orogen"]
    halo_cells = {"Continental Arc": 3, "Island Arc": 3, "Orogen": 3,
                  "Continental Margin": 3, "Extended Crust": 4,
                  "Basin": 3, "Platform": 6, "Shield": 8}
    for prov in paint_priority:
        sub = dec[dec["prov"] == prov]
        if sub.empty:
            continue
        h = halo_cells.get(prov, 3)
        for _, r in sub.iterrows():
            iy, ix = int(r["iy"]), int(r["ix"])
            i0, i1 = max(0, iy-h), min(nlat, iy+h+1)
            j0, j1 = max(0, ix-h), min(nlon, ix+h+1)
            block_cont = cont_mask[i0:i1, j0:j1]
            P_grid[i0:i1, j0:j1] = np.where(block_cont, PROV_INDEX[prov],
                                            P_grid[i0:i1, j0:j1])
    return P_grid


if __name__ == "__main__":
    # quick sanity check
    import time
    lat = np.linspace(-90, 90, 181)
    lon = np.linspace(-180, 180, 361)
    from sw_io import load_grid
    for t in [0, 50, 100, 250, 500]:
        try:
            M, _, _ = load_grid(t)
        except FileNotFoundError as e:
            print(f"t={t}: skip ({e})")
            continue
        t0 = time.time()
        mask = cob_mask(t, lat, lon, M)
        sz = subduction_zones(t)
        d = distance_to_subduction(t, lat, lon)
        print(f"t={t} Ma: cont cells={mask.sum()}, "
              f"SZ segs={len(sz)} verts={sum(len(s) for s in sz)}, "
              f"d_sz_min={d.min():.0f} km, {time.time()-t0:.1f}s")
