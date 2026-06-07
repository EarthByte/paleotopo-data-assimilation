"""
=============================================================================
hasterok_provinces.py — Alternative tectonic-province classifier
=============================================================================

Drop-in replacement for the geomorphic-surrogate `province_grid()` in
`plate_model_utils_scotese.py`, using the global geological-province
shapefile of:

    Hasterok, D., Halpin, J.A., Collins, A.S., Hand, M., Kreemer, C.,
    Gard, M.G. & Glorie, S., 2022.  New maps of global geological
    provinces and tectonic plates.  Earth-Science Reviews, 231, 104069.

Shapefile expected at:
    data/Hasterok_plates_provinces/shp/global_gprv.shp

The Hasterok dataset gives a finer-grained classification (15 province
types, 914 polygons) than the existing 9-category Scotese scheme:

    volcanic arc, orogenic belt, accretionary complex, shield,
    oceanic crust, craton, passive margin, wide rift, back-arc basin,
    oceanic back-arc basin, narrow rift, ophiolite complex,
    magmatic province, foredeep basin, basin
    (plus an "Other" fall-through for unclassified cells)

USAGE
    from hasterok_provinces import (
        HASTEROK_PROV_LIST, HASTEROK_PROV_INDEX, HASTEROK_POOLS,
        province_grid_hasterok, classify_samples_hasterok,
    )

    P_g = province_grid_hasterok(t, lat, lon, cont_mask, dec)
    df["Tecto_Prov"] = classify_samples_hasterok(df["Lat"], df["Lon"])

ONE-TIME SETUP (auto-cached on first call)
    - Each polygon centroid is assigned a S&W plate ID by point-in-polygon
      against the Scotese 2008 present-day continental polygons.
    - Per-sample classification is cached to .npy alongside the polygon
      plate-ID cache.

DEPENDENCIES
    geopandas, shapely, rasterio (for fast polygon rasterisation),
    pygplates (for reconstruction), numpy.
=============================================================================
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np

import geopandas as gpd
import pygplates as pg
import rasterio.features
from shapely.geometry import Polygon, MultiPolygon, mapping, shape
from shapely.affinity import translate
from rasterio.transform import from_bounds

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from paths_scotese import PLATE_DIR

# ---------------------------------------------------------------------------
# Province list + integer codes
# ---------------------------------------------------------------------------
HASTEROK_PROV_LIST = [
    "volcanic arc",
    "orogenic belt",
    "accretionary complex",
    "shield",
    "oceanic crust",
    "craton",
    "passive margin",
    "wide rift",
    "back-arc basin",
    "oceanic back-arc basin",
    "narrow rift",
    "ophiolite complex",
    "magmatic province",
    "foredeep basin",
    "basin",
    "Other",
]
HASTEROK_PROV_INDEX = {p: i for i, p in enumerate(HASTEROK_PROV_LIST)}

# Pooling table for the smooth-shrinkage refinement.  Group related
# tectonic settings so a province with few samples can borrow CDF
# support from its neighbours.
HASTEROK_POOLS = {
    "volcanic arc":           ["volcanic arc", "accretionary complex", "orogenic belt"],
    "orogenic belt":          ["orogenic belt", "accretionary complex", "volcanic arc",
                               "magmatic province"],
    "accretionary complex":   ["accretionary complex", "orogenic belt", "volcanic arc"],
    "shield":                 ["shield", "craton"],
    "craton":                 ["craton", "shield"],
    "passive margin":         ["passive margin", "ophiolite complex",
                               "wide rift", "narrow rift"],
    "wide rift":              ["wide rift", "narrow rift", "passive margin"],
    "narrow rift":            ["narrow rift", "wide rift"],
    "back-arc basin":         ["back-arc basin", "foredeep basin", "basin"],
    "oceanic back-arc basin": ["oceanic back-arc basin", "back-arc basin"],
    "ophiolite complex":      ["ophiolite complex", "passive margin",
                               "accretionary complex"],
    "magmatic province":      ["magmatic province", "orogenic belt"],
    "foredeep basin":         ["foredeep basin", "basin", "back-arc basin"],
    "basin":                  ["basin", "foredeep basin"],
    "oceanic crust":          ["oceanic crust"],
    "Other":                  ["Other"],
}

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SHP_PATH = (Path(__file__).resolve().parent.parent /
            "data" / "Hasterok_plates_provinces" / "shp" / "global_gprv.shp")
PLATE_ID_CACHE = SHP_PATH.parent / "hasterok_polygon_plate_ids.npy"
CONTINENTAL_POLYGONS_FILE = (
    PLATE_DIR / "Scotese_Wright_PresentDay_ContinentalPolygons.gpml"
)

# ---------------------------------------------------------------------------
# Polygon loading + plate-ID assignment (one-time, cached)
# ---------------------------------------------------------------------------
_GDF = None              # geopandas GeoDataFrame, cached in memory
_ROT_MODEL = None        # pyGPlates RotationModel, cached in memory


def _rotation_model():
    global _ROT_MODEL
    if _ROT_MODEL is None:
        from paths_scotese import ROT_FILE
        _ROT_MODEL = pg.RotationModel(str(ROT_FILE))
    return _ROT_MODEL


def load_polygons():
    """Load the Hasterok shapefile once, cached in memory.  Assigns each
    polygon a S&W plate ID at present-day, caching the result to .npy."""
    global _GDF
    if _GDF is not None:
        return _GDF
    if not SHP_PATH.exists():
        raise FileNotFoundError(
            f"Hasterok shapefile not found at {SHP_PATH}.  "
            "Download from Hasterok et al. (2022, ESR) supplementary."
        )
    gdf = gpd.read_file(SHP_PATH)
    if "prov_type" not in gdf.columns:
        raise ValueError(
            f"Shapefile has no 'prov_type' column.  Columns: {list(gdf.columns)}"
        )
    # Fill any missing prov_type with "Other"
    gdf["prov_type"] = gdf["prov_type"].fillna("Other").astype(str)
    # Quantize to known categories — anything outside the list goes to "Other"
    known = set(HASTEROK_PROV_LIST)
    gdf["prov_type"] = gdf["prov_type"].where(gdf["prov_type"].isin(known), "Other")

    # Plate-ID assignment
    if PLATE_ID_CACHE.exists():
        pids = np.load(PLATE_ID_CACHE)
        if len(pids) == len(gdf):
            gdf["sw_plate_id"] = pids
            _GDF = gdf
            return gdf
    print(f"  [hasterok] assigning S&W plate IDs to {len(gdf)} polygons …")
    if not CONTINENTAL_POLYGONS_FILE.exists():
        print(f"  [hasterok] continental polygons file not found at "
              f"{CONTINENTAL_POLYGONS_FILE}; defaulting all plate IDs to 0")
        gdf["sw_plate_id"] = 0
    else:
        partitioner = pg.PlatePartitioner(str(CONTINENTAL_POLYGONS_FILE),
                                          _rotation_model())
        plate_ids = []
        for geom in gdf.geometry:
            # Use polygon centroid as the assignment point.
            c = geom.representative_point()
            point = pg.PointOnSphere(c.y, c.x)
            partition = partitioner.partition_point(point)
            if partition is None:
                plate_ids.append(0)
            else:
                pid = partition.get_feature().get_reconstruction_plate_id()
                plate_ids.append(pid if pid is not None else 0)
        gdf["sw_plate_id"] = np.array(plate_ids, dtype=np.int32)
        np.save(PLATE_ID_CACHE, gdf["sw_plate_id"].to_numpy())
        print(f"  [hasterok] cached plate IDs to {PLATE_ID_CACHE}")
    _GDF = gdf
    return gdf


# ---------------------------------------------------------------------------
# Per-age reconstruction + rasterisation
# ---------------------------------------------------------------------------
def _reconstruct_one(geom, plate_id, t, rot_model):
    """Reconstruct a single shapely Polygon/MultiPolygon to age t via
    plate_id under rot_model.  Returns a shapely geometry (or None if
    the rotation isn't defined or the polygon is empty)."""
    if plate_id is None or plate_id == 0 or t == 0:
        return geom
    try:
        rot = rot_model.get_rotation(float(t), int(plate_id),
                                     anchor_plate_id=0)
    except Exception:
        return None
    if rot is None:
        return geom

    def _rot_ring(coords):
        out = []
        for x, y in coords:
            pp = pg.PointOnSphere(y, x)
            rp = rot * pp
            ll = rp.to_lat_lon()
            out.append((ll[1], ll[0]))
        return out

    def _rot_polygon(p: Polygon) -> Polygon:
        ext = _rot_ring(list(p.exterior.coords))
        ints = [_rot_ring(list(r.coords)) for r in p.interiors]
        return Polygon(ext, ints)

    if isinstance(geom, Polygon):
        return _rot_polygon(geom)
    if isinstance(geom, MultiPolygon):
        return MultiPolygon([_rot_polygon(p) for p in geom.geoms])
    return None


def province_grid_hasterok(t, lat, lon, cont_mask):
    """Build a (nlat, nlon) int8 raster of Hasterok province codes at
    age t.  Continental cells outside any reconstructed polygon get the
    "Other" index; non-continental cells get -1."""
    gdf = load_polygons()
    rot = _rotation_model()

    nlat, nlon = len(lat), len(lon)
    P_grid = np.full((nlat, nlon), HASTEROK_PROV_INDEX["Other"], dtype=np.int8)
    P_grid[~cont_mask] = -1

    # rasterio transform: lon increases west→east, lat increases south→north.
    # Hasterok bounds: lon in [-180, 180], lat in [-90, 90].
    transform = from_bounds(lon[0] - 0.5, lat[0] - 0.5,
                            lon[-1] + 0.5, lat[-1] + 0.5,
                            nlon, nlat)

    # Build (geom, code) pairs.  Lower-priority types are rasterised first
    # so higher-priority types over-paint them.  Priority order (from
    # broadest to most localised) is the reverse of the listed precedence:
    # cratons / shields paint LAST so they "win" in the interior.  Actually
    # for province classification we want the most specific tectonic type
    # to win, so we paint cratons/shields FIRST (broad cratonic background)
    # and then over-paint with arcs / orogens / rifts / basins.
    PRIORITY = [
        "oceanic crust", "craton", "shield", "passive margin",
        "magmatic province", "basin", "foredeep basin",
        "wide rift", "narrow rift", "ophiolite complex",
        "back-arc basin", "oceanic back-arc basin",
        "accretionary complex", "orogenic belt", "volcanic arc",
        "Other",
    ]

    for prov in PRIORITY:
        sub = gdf[gdf["prov_type"] == prov]
        if sub.empty:
            continue
        shapes = []
        code = HASTEROK_PROV_INDEX[prov]
        for _, row in sub.iterrows():
            reco = _reconstruct_one(row.geometry,
                                    int(row["sw_plate_id"]),
                                    float(t), rot)
            if reco is None or reco.is_empty:
                continue
            shapes.append((mapping(reco), code))
        if not shapes:
            continue
        try:
            mask = rasterio.features.rasterize(
                shapes,
                out_shape=(nlat, nlon),
                transform=transform,
                fill=-1,
                dtype="int16",
                all_touched=False,
            )
        except Exception as e:
            print(f"  [hasterok] rasterise failed for '{prov}' at {t} Ma: {e}")
            continue
        # Only over-write continental cells; preserve -1 for ocean.
        update = (mask >= 0) & cont_mask
        P_grid[update] = mask[update].astype(np.int8)

    return P_grid


# ---------------------------------------------------------------------------
# Sample classification at present-day
# ---------------------------------------------------------------------------
_SAMPLE_PROV_CACHE = SHP_PATH.parent / "hasterok_sample_prov_labels.npz"


def classify_samples_hasterok(lat_arr, lon_arr, cache_key=None):
    """Return a numpy array of Hasterok province-type strings for each
    (lat, lon) sample.  Uses a small spatial index for speed.  Cache key
    can be the length of the input arrays (used to validate the cache)."""
    n = len(lat_arr)
    if _SAMPLE_PROV_CACHE.exists():
        d = np.load(_SAMPLE_PROV_CACHE, allow_pickle=True)
        if int(d["n"]) == n:
            return d["labels"]
    print(f"  [hasterok] classifying {n} samples by point-in-polygon …")
    gdf = load_polygons()
    from shapely.geometry import Point
    from shapely.strtree import STRtree
    polys = list(gdf.geometry)
    types = list(gdf["prov_type"])
    tree = STRtree(polys)
    out = np.full(n, "Other", dtype=object)
    for i, (la, lo) in enumerate(zip(lat_arr, lon_arr)):
        if not (np.isfinite(la) and np.isfinite(lo)):
            continue
        pt = Point(lo, la)
        # Query candidates from the spatial index; STRtree.query returns
        # indices in shapely 2.x, geometries in 1.x — handle both.
        candidates = tree.query(pt)
        for ci in candidates:
            poly = polys[ci] if isinstance(ci, (int, np.integer)) else ci
            if poly.contains(pt):
                t_idx = types[ci] if isinstance(ci, (int, np.integer)) \
                        else types[polys.index(ci)]
                out[i] = t_idx
                break
    np.savez(_SAMPLE_PROV_CACHE, n=n, labels=out)
    print(f"  [hasterok] cached sample labels to {_SAMPLE_PROV_CACHE}")
    return out


if __name__ == "__main__":
    # Sanity smoke test
    gdf = load_polygons()
    print(f"loaded {len(gdf)} polygons, {gdf['prov_type'].nunique()} unique types")
    print("plate-ID distribution (top 10):")
    print(gdf["sw_plate_id"].value_counts().head(10))
