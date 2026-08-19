#!/usr/bin/env python3
"""shirmard_cratons.py — seismically defined craton mask for sample selection.

Cratonic surface elevation is not set by crustal thickness.  Mooney and Vidale
(2003, Lithos 71, 185-193) show that cratons sit 300 to 1500 m below the
elevation predicted by pure crustal isostasy, because the compensation is
carried by sub-crustal density -- roots that are colder, or less depleted, than
thermal and xenolith models assume.  A chemical-mohometry crustal thickness
therefore has no skill in predicting cratonic elevation, and an Airy conversion
applied there is biased high, which would add spurious relief and spurious
erosion to the corrected paleotopography.

Craton outlines here come from

    Shirmard, H., et al. (2025). Craton boundary detection from full-waveform
    tomography model reveals links to critical metal deposits.
    Geoscience Frontiers 16, 102176.  doi:10.1016/j.gsf.2025.102176
    Data: doi:10.5281/zenodo.15003746

specifically ``Craton_Boundaries_Inferred.shp``, 20 named craton polygons
totalling 5.00e7 km2.  The companion ``Craton_Boundaries.shp`` holds the 546
detected boundary line segments and is not used: it is not a closed mask.

These outlines are preferred over the craton and shield classes of the Hasterok
et al. (2022) province compilation because they are defined by lithospheric
structure in a full-waveform tomography model -- thick, fast, seismically
distinct roots -- which is the same physical property that makes cratonic
elevation non-isostatic in the crustal sense.  The Hasterok classes are
defined by surface geology and last orogenic event.  The two agree on much
less than their similar total areas suggest: 3.04e7 km2 in common against
5.00e7 (Shirmard) and 5.38e7 (Hasterok craton plus shield), a Jaccard index
of 0.42.

Samples are tested at their PRESENT-DAY coordinates, matching how
``assimilate_scotese.prepare_samples`` assigns tectonic provinces.  No plate
reconstruction of these polygons is required, and none should be added: the
mask selects which observations enter the correction, and never acts on the
correction field itself.

A polygon mask is also the only workable route here.  ``PROVINCE_CLASSIFIER``
defaults to ``"geomorphic"``, whose nine classes include ``Shield`` but no
craton class at all, so excluding cratons by province name would depend on
which classifier happened to be active and would silently miss most cratonic
area under the default.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths_scotese import PROJECT_ROOT  # noqa: E402

# PROJECT_ROOT honours PALEOTOPO_PROJECT_ROOT and points at the repository
# root, whose data/ directory holds every other input (LithoRef18, the plate
# model, the Scotese DEMs).  Resolving the path from __file__ instead would
# land in public_repo_scotese/data/, which is a different directory and does
# not hold these inputs.
SHP_PATH = (PROJECT_ROOT / "data" / "Shirmard_cratons" /
            "Craton_Boundaries_Inferred.shp")

_POLY_CACHE = None
_TREE_CACHE = None


def load_craton_polygons():
    """Load and cache the Shirmard craton polygons as a GeoDataFrame."""
    global _POLY_CACHE
    if _POLY_CACHE is None:
        import geopandas as gpd
        if not SHP_PATH.exists():
            raise FileNotFoundError(
                f"Shirmard craton polygons not found at {SHP_PATH}. "
                "Copy Craton_Boundaries_Inferred.{shp,shx,dbf,prj} from the "
                "Shirmard et al. (2025) Zenodo archive into data/Shirmard_cratons/."
            )
        g = gpd.read_file(SHP_PATH)
        if g.crs is not None and g.crs.to_epsg() != 4326:
            g = g.to_crs(epsg=4326)
        _POLY_CACHE = g
    return _POLY_CACHE


def _tree():
    global _TREE_CACHE
    if _TREE_CACHE is None:
        from shapely.strtree import STRtree
        g = load_craton_polygons()
        _TREE_CACHE = (STRtree(list(g.geometry)), list(g.geometry))
    return _TREE_CACHE


def in_craton(lat_arr, lon_arr) -> np.ndarray:
    """Boolean mask: True where the (lat, lon) point falls inside a craton.

    Longitudes are accepted in either [-180, 180] or [0, 360] and wrapped to
    the polygons' [-180, 180] convention.  Non-finite coordinates return False.
    """
    from shapely.geometry import Point
    lat = np.asarray(lat_arr, dtype=float)
    lon = np.asarray(lon_arr, dtype=float)
    if lat.shape != lon.shape:
        raise ValueError("lat and lon must have the same shape")
    lon = ((lon + 180.0) % 360.0) - 180.0
    out = np.zeros(lat.shape, dtype=bool)
    tree, polys = _tree()
    finite = np.isfinite(lat) & np.isfinite(lon)
    for i in np.flatnonzero(finite.ravel()):
        idx = np.unravel_index(i, lat.shape)
        pt = Point(float(lon[idx]), float(lat[idx]))
        for cand in tree.query(pt):
            # shapely 2.x returns indices, 1.x returns geometries
            geom = polys[cand] if isinstance(cand, (int, np.integer)) else cand
            if geom.contains(pt):
                out[idx] = True
                break
    return out


def craton_area_km2() -> float:
    """Total area of the craton mask, for reporting."""
    g = load_craton_polygons()
    return float(g.geometry.to_crs("+proj=cea").area.sum() / 1e6)


def _self_test() -> int:
    g = load_craton_polygons()
    print(f"[shirmard] {len(g)} craton polygons, {craton_area_km2():.3e} km2")
    # Points chosen well inside and well outside named cratons.
    inside = [(-26.0, 120.0, "Yilgarn"), (-27.0, 26.0, "Kaapvaal"),
              (60.0, -100.0, "North American"), (62.0, 105.0, "Siberian")]
    outside = [(-33.9, 151.2, "Sydney"), (35.7, 139.7, "Tokyo"),
               (0.0, -140.0, "central Pacific"), (46.0, 8.0, "Alps")]
    ok = True
    for la, lo, nm in inside:
        r = bool(in_craton([la], [lo])[0])
        ok &= r
        print(f"  {'PASS' if r else 'FAIL'} inside  {nm:<16} ({la}, {lo}) -> {r}")
    for la, lo, nm in outside:
        r = bool(in_craton([la], [lo])[0])
        ok &= (not r)
        print(f"  {'PASS' if not r else 'FAIL'} outside {nm:<16} ({la}, {lo}) -> {r}")
    # Longitude wrapping must agree
    a = in_craton([62.0], [105.0])[0]
    b = in_craton([62.0], [465.0])[0]
    ok &= (a == b)
    print(f"  {'PASS' if a == b else 'FAIL'} longitude wrap 105 == 465")
    print("ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(_self_test())
