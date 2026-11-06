#!/usr/bin/env python3
"""esr_revision_diagnostics.py — per-age quantities of the corrected grids, computed once: the support layer, the
declustered observations with their in-sample residuals, and the area above elevation thresholds.

For every age the workflow's own assimilate_one runs in memory with the product's settings (CONTROL: it must
reproduce the grid in CORRECTED_DIR to 1 m, or the script stops). Written to OUTPUT_DIR/esr_revision, one file each:

  support_km.nc        support_km(age, lat, lon): great-circle distance from every cell to the nearest declustered
                       observation of that age (km; NaN where the age has none), computed by the workflow's own
                       nearest_sample_distance_grid, so it equals the support_km variable of the grids. The residual
                       kernel is whole within 250 km and tapered to zero at 500 km; the province rescaling is zero
                       beyond 400 km.
  declustered.csv      every declustered observation: age, reconstructed position, cell, province, proxy elevation z,
                       sigma, the input and corrected elevations at its cell.
  per_age.csv          per age: number of declustered observations, the area fraction of continental cells above
                       2 and 3 km in the input and corrected grids (cos-latitude weights), land p99 before and after,
                       continental dz RMS, the share of continental area farther than 500 km from any observation.

    cd scripts
    python3 esr_revision_diagnostics.py                 # all ages, ~2 min
    python3 esr_revision_diagnostics.py --ages 100,500  # quick check
"""
from __future__ import annotations

import argparse
import contextlib
import io
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from paths_scotese import CORRECTED_DIR as PRODUCT, OUTPUT_DIR  # noqa: E402

OUT = OUTPUT_DIR / "esr_revision"
R_KM = 6371.0


def xyz(lat, lon):
    la, lo = np.radians(lat), np.radians(lon)
    return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], -1)


def support_grid(lat, lon, slat, slon):
    """Distance (km) from every cell centre to the nearest of the points (slat, slon); NaN without points."""
    from scipy.spatial import cKDTree
    LA, LO = np.meshgrid(lat, lon, indexing="ij")
    if len(slat) == 0:
        return np.full(LA.shape, np.nan, np.float32)
    d, _ = cKDTree(xyz(np.asarray(slat), np.asarray(slon))).query(xyz(LA.ravel(), LO.ravel()))
    return (2.0 * R_KM * np.arcsin(np.clip(d / 2.0, 0.0, 1.0))).reshape(LA.shape).astype(np.float32)


def area_frac(field, cont, lat, thr):
    w = np.cos(np.radians(lat))[:, None] * np.ones((1, field.shape[1]))
    ca = (cont * w).sum()
    return float(((cont & (field > thr)) * w).sum() / ca) if ca > 0 else np.nan


def _self_test():
    lat = np.array([0.0, 10.0]); lon = np.array([0.0, 10.0])
    g = support_grid(lat, lon, [0.0], [0.0])
    assert abs(g[0, 0]) < 1e-6 and abs(g[1, 0] - 1111.95) < 1.0, g
    assert np.isnan(support_grid(lat, lon, [], [])).all()
    cont = np.array([[True, True], [False, True]]); f = np.array([[3500.0, 100.0], [9999.0, 2500.0]])
    assert 0.3 < area_frac(f, cont, lat, 3000.0) < 0.35
    print("[esr-diag] self-test passed")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--ages", default="all")
    ap.add_argument("--out", type=Path, default=OUT)
    a = ap.parse_args()
    _self_test()
    if a.self_test:
        return 0
    import pandas as pd
    import netCDF4 as nc
    import assimilate_scotese as A
    from sw_io import available_ages
    ages = [int(round(t)) for t in available_ages()] if a.ages == "all" else [int(x) for x in a.ages.split(",")]
    a.out.mkdir(parents=True, exist_ok=True)
    rows, decs, sup = [], [], []
    lat = lon = None
    for t in ages:
        with contextlib.redirect_stdout(io.StringIO()):
            s, M, Mf, P, cont, dec = A.assimilate_one(float(t), save_nc=False)
        _, lat, lon = A.sw_load_grid(float(t))
        with nc.Dataset(PRODUCT / f"{t}Ma_corrected_SW.nc") as d:
            Mr = np.ma.filled(d["M_corrected"][:], np.nan).astype(float)
            cmask = np.ma.filled(d["continent_mask"][:], 0) > 0
        dmax = float(np.nanmax(np.abs(Mr - Mf)))
        if not dmax <= 1.0:
            raise SystemExit(f"[esr-diag] CONTROL FAILED at {t} Ma: {dmax:.2f} m from {PRODUCT.name}")
        if len(dec):
            sl = np.asarray(lat)[dec["iy"].astype(int).values]; so = np.asarray(lon)[dec["ix"].astype(int).values]
        else:
            sl = so = np.array([])
        g = A.nearest_sample_distance_grid(lat, lon, dec).astype(np.float32)   # assimilate_one runs with save_nc=False
        g[~np.isfinite(g)] = np.nan
        sup.append(g)
        if len(dec):
            iy, ix = dec["iy"].astype(int).values, dec["ix"].astype(int).values
            decs.append(pd.DataFrame(dict(age=t, rlat=dec["rlat"].values, rlon=dec["rlon"].values, cell_lat=sl, cell_lon=so,
                                          prov=dec["prov"].values, z=dec["z"].values, sigma=dec["sigma"].values,
                                          prior=M[iy, ix], corrected=Mf[iy, ix])))
        land0, land1 = cmask & (M > 0), cmask & (Mf > 0)
        w = np.cos(np.radians(lat))[:, None] * np.ones((1, len(lon)))
        rows.append(dict(age=t, n_decluster=int(len(dec)),
                         frac_gt2_input=area_frac(M, cmask, lat, 2000.0), frac_gt3_input=area_frac(M, cmask, lat, 3000.0),
                         frac_gt2_corrected=area_frac(Mf, cmask, lat, 2000.0), frac_gt3_corrected=area_frac(Mf, cmask, lat, 3000.0),
                         p99_input=float(np.percentile(M[land0], 99)) if land0.any() else np.nan,
                         p99_corrected=float(np.percentile(Mf[land1], 99)) if land1.any() else np.nan,
                         dz_rms=float(np.sqrt(np.nanmean((Mf - M)[cmask] ** 2))),
                         frac_cont_beyond_500km=float(((cmask & ~(g <= 500.0)) * w).sum() / max((cmask * w).sum(), 1e-12))))
        print(f"[esr-diag] {t:3d} Ma: {len(dec):3d} declustered; control {dmax:.3f} m; > 3 km {100 * rows[-1]['frac_gt3_corrected']:.2f} %; "
              f"continent beyond 500 km {100 * rows[-1]['frac_cont_beyond_500km']:.0f} %", flush=True)
    pd.DataFrame(rows).to_csv(a.out / "per_age.csv.partial", index=False, float_format="%.6g")
    os.replace(a.out / "per_age.csv.partial", a.out / "per_age.csv")
    pd.concat(decs).to_csv(a.out / "declustered.csv.partial", index=False, float_format="%.4f")
    os.replace(a.out / "declustered.csv.partial", a.out / "declustered.csv")
    f = a.out / "support_km.nc.partial"
    with nc.Dataset(f, "w") as d:
        d.createDimension("age", len(ages)); d.createDimension("lat", len(lat)); d.createDimension("lon", len(lon))
        d.createVariable("age", "f4", ("age",))[:] = ages
        d.createVariable("lat", "f4", ("lat",))[:] = lat
        d.createVariable("lon", "f4", ("lon",))[:] = lon
        v = d.createVariable("support_km", "f4", ("age", "lat", "lon"), zlib=True, complevel=4)
        v[:] = np.stack(sup); v.units = "km"
        v.long_name = "great-circle distance to the nearest declustered geochemical observation of that age"
        d.title = "Support layer of the corrected Scotese-Wright PaleoDEMs"
        d.note = ("residual kernel whole within 250 km of the nearest observation, tapered to zero at 500 km; "
                  "province rescaling zero beyond 400 km; beyond ~500 km the corrected grid is the input grid")
    os.replace(f, a.out / "support_km.nc")
    print(f"[esr-diag] wrote {a.out}/per_age.csv, declustered.csv, support_km.nc")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
