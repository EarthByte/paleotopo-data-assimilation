#!/usr/bin/env python3
"""east_asia_ablation.py — which settings of the Scotese-Wright correction shape eastern Asia, and which configuration
the independent paleoelevation records there support.

THE QUESTION. Against the first version of the correction ("submitted" below), the product of record (1) puts some
East Asian basins lower, (2) breaks the > 3 km terrain into separate bodies and (3) at 120 Ma makes eastern North China
higher than the west. The continental mask on the sphere, the craton exclusion and the young-age taper are common to
every configuration; the ablation asks what the other settings do: the reach limit (residual kernel whole to 250 km,
zero at 500 km), the support limit (province rescaling to 400 km), the final smoothing pass, the cross-boundary
smoothing and the removal of samples on the input map's sea.

THE ABLATION. assimilate_scotese.assimilate_one runs in memory (nothing is written to a product directory) for each age
and configuration. v110 is the module's own settings; every other configuration reverts one or more changes:
  v110              the product of record (CONTROL: must reproduce the grids in CORRECTED_DIR to 1 m at every age)
  rev_support       LOCAL_SUPPORT_KM 1,000 (submitted)
  rev_reach         RESID_REACH_KM off (submitted)
  rev_reach_support both limits reverted
  rev_final         the final 4.5-cell pass on dz back on (submitted)
  rev_cross         cross-boundary smoothing 3.0 cells (submitted)
  rev_sea           samples on the input map's sea kept (submitted)
  rev_craton        craton samples kept (submitted)
  submitted_fixed   all of the above reverted except the craton exclusion (the spherical mask is in every
                    configuration)
  submitted         submitted_fixed with the craton samples kept (the first version of the method, on the spherical
                    mask)

MEASURED per configuration and age:
  records   every independent record of data/independent_records/east_asia_paleoelevation_200-66Ma.csv (columns lo_m /
            hi_m give the interval a record admits, sea-level datums -250..+250 m; `use` strong / weak / tibet / added,
            fixed before any score was computed) at every grid age inside its age interval (or the nearest grid age
            within 5 Myr of its midpoint): the site reconstructed with the workflow's own reconstructor
            (ScoteseSampleReconstructor), the median elevation of the 3 x 3 cells around it, the misfit
            max(lo - z, 0, z - hi), and the site's distance to the nearest declustered training sample. Per record the
            misfit is averaged over its ages, so every record counts once.
  contrast  median elevation of continental cells East minus West (present-day boxes, reconstructed; defaults West
            34-41N 104-112E, the Ordos and Taihang; East 30-41N 114-123E, eastern North China, Hefei, Shandong).
  high      cells above 3,000 m within the East Asia box (20-55N, 95-145E present day, reconstructed): connected bodies
            of at least 3 cells (8-connected), the largest body's and the total area (km2).
  basins    median elevation of the 3 x 3 cells at Songliao, Erlian, Hailar, Ordos, Sichuan, Hefei, Jiaolai, Jixi.

DECISION RULE (fixed before any score was computed). R1, on the STRONG records (use == strong; Tibet and the weak
records are reported, not ruled): the reach and support limits stay if the mean misfit of v110 is no more than 100 m
above the better of rev_reach_support and submitted_fixed. If either fits the strong records better by more than
100 m, the far-field rule is revised globally (never for eastern Asia alone). R1 is also reported on the strong
records farther than 150 km from the nearest training sample, where the two limits act. The contrast, the high bodies
and the basins attribute each of the three observations to a single revert; they are reported, not ruled.

    cd scripts
    python3 east_asia_ablation.py --ages 120 --jobs 1                 # pilot: one age, every configuration
    python3 east_asia_ablation.py --jobs 6                            # 200-65 Ma, every 5 Myr
"""
from __future__ import annotations

import argparse
import contextlib
import io
import json
import os
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from paths_scotese import PROJECT_ROOT as ROOT, CORRECTED_DIR as CONTROL_DIR, OUTPUT_DIR  # noqa: E402

RECORDS = ROOT / "data" / "independent_records" / "east_asia_paleoelevation_200-66Ma.csv"
SUB = {"LOCAL_SUPPORT_KM": 1000.0, "RESID_REACH_KM": (0.0, 0.0), "DELTA_SMOOTH_SIGMA_FINAL": 4.5,
       "CROSS_BOUNDARY_SIGMA": 3.0, "DROP_SEA_SITED_SAMPLES": False}
CONFIGS = {
    "v110": {},
    "rev_support": {"LOCAL_SUPPORT_KM": 1000.0},
    "rev_reach": {"RESID_REACH_KM": (0.0, 0.0)},
    "rev_reach_support": {"LOCAL_SUPPORT_KM": 1000.0, "RESID_REACH_KM": (0.0, 0.0)},
    "rev_final": {"DELTA_SMOOTH_SIGMA_FINAL": 4.5},
    "rev_cross": {"CROSS_BOUNDARY_SIGMA": 3.0},
    "rev_sea": {"DROP_SEA_SITED_SAMPLES": False},
    "rev_craton": {"DROP_CRATON_SAMPLES": False},
    "submitted_fixed": dict(SUB),
    "submitted": {**SUB, "DROP_CRATON_SAMPLES": False},
}
EA_BOX = (20.0, 55.0, 95.0, 145.0)
WEST_BOX = (34.0, 41.0, 104.0, 112.0)
EAST_BOX = (30.0, 41.0, 114.0, 123.0)
BASINS = {"Songliao": (45.5, 124.8), "Erlian": (43.5, 112.5), "Hailar": (49.0, 119.5), "Ordos": (37.5, 108.5),
          "Sichuan": (30.5, 105.5), "Hefei": (32.0, 117.0), "Jiaolai": (36.9, 120.7), "Jixi": (45.3, 131.0)}
HIGH_M, MIN_BODY = 3000.0, 3
TOL_M, FAR_KM = 100.0, 150.0
R_KM = 6371.0


def misfit(z, lo, hi):
    """Distance of z outside [lo, hi] (0 inside). Pure."""
    return float(max(lo - z, 0.0, z - hi))


def record_ages(amin, amax, grid_ages):
    """Grid ages inside [amin, amax]; else the nearest within 5 Myr of the midpoint; else none. Pure."""
    a = [t for t in grid_ages if amin - 1e-9 <= t <= amax + 1e-9]
    if a:
        return a
    mid = 0.5 * (amin + amax)
    near = min(grid_ages, key=lambda t: abs(t - mid)) if grid_ages else None
    return [near] if near is not None and abs(near - mid) <= 5.0 else []


def bodies(mask, min_cells=MIN_BODY):
    """8-connected bodies of a boolean grid with at least min_cells cells: list of index arrays. Pure."""
    from scipy import ndimage
    lab, n = ndimage.label(mask, structure=np.ones((3, 3), bool))
    out = []
    for k in range(1, n + 1):
        idx = np.where(lab == k)
        if len(idx[0]) >= min_cells:
            out.append(idx)
    return out


def _self_test():
    assert misfit(500, 0, 1000) == 0 and misfit(-300, -250, 250) == 50 and misfit(3000, 1500, np.inf) == 0
    assert misfit(1000, -np.inf, 600) == 400
    assert record_ages(113.2, 121.4, [110, 115, 120, 125]) == [115, 120]
    assert record_ages(132.17, 132.49, [125, 130, 135]) == [130]
    assert record_ages(83.72, 84.72, [70, 75]) == []
    m = np.zeros((6, 6), bool); m[0, 0:3] = True; m[4, 4] = True; m[3, 3] = True; m[5, 5] = True
    b = bodies(m)
    assert len(b) == 2 and sorted(len(x[0]) for x in b) == [3, 3], [len(x[0]) for x in b]   # diagonal run joins
    print("[ea] self-test passed")


def gc_km(la1, lo1, la2, lo2):
    la1, lo1, la2, lo2 = map(np.radians, (la1, lo1, la2, lo2))
    a = np.sin((la2 - la1) / 2) ** 2 + np.cos(la1) * np.cos(la2) * np.sin((lo2 - lo1) / 2) ** 2
    return 2 * R_KM * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def box_points(b, step=0.5):
    la = np.arange(b[0], b[1] + 1e-9, step); lo = np.arange(b[2], b[3] + 1e-9, step)
    LA, LO = np.meshgrid(la, lo, indexing="ij")
    return LA.ravel(), LO.ravel()


def worker(args):
    ages, names, recs, west_box, east_box = args
    os.environ["PALEO_CORRECTED_DIR"] = tempfile.mkdtemp(prefix="ea_ablation_")
    os.environ["PALEO_OUTPUT_DIR"] = os.environ["PALEO_CORRECTED_DIR"]
    import importlib
    import paths_scotese
    importlib.reload(paths_scotese)                     # assimilate_scotese, imported next, writes nowhere near the product
    import pandas as pd
    import netCDF4 as nc
    import assimilate_scotese as A
    from sw_io import nearest_cell_index
    defaults = {k: getattr(A, k) for v in CONFIGS.values() for k in v}
    sr = A.get_reconstructor()

    def frame(lat, lon):                                # present-day points, plate ids assigned once
        df = pd.DataFrame({"Lat": np.asarray(lat, float), "Lon": np.asarray(lon, float)})
        df["sw_plate_id"] = sr.assign_plate_ids(df)
        return df

    F = {"EA": frame(*box_points(EA_BOX)), "W": frame(*box_points(west_box)), "E": frame(*box_points(east_box)),
         "sites": frame(recs["lat"].values, recs["lon"].values),
         "basins": frame([v[0] for v in BASINS.values()], [v[1] for v in BASINS.values()])}

    def recon(key, t):
        return sr.reconstruct(F[key], float(t))
    rec_rows, age_rows, basin_rows = [], [], []
    for t in ages:
        M0, lat, lon = A.sw_load_grid(float(t))
        cell_km2 = (111.195 ** 2) * np.cos(np.radians(lat))[:, None] * np.ones((1, len(lon)))

        def cells(la, lo):
            ok = np.isfinite(la) & np.isfinite(lo)
            return nearest_cell_index(lat, la[ok]), nearest_cell_index(lon, lo[ok])

        R = {k: recon(k, t) for k in ("EA", "W", "E")}
        site_lat, site_lon = recon("sites", t)
        b_lat, b_lon = recon("basins", t)
        for name in names:
            for k, v in defaults.items():
                setattr(A, k, v)
            for k, v in CONFIGS[name].items():
                setattr(A, k, v)
            with contextlib.redirect_stdout(io.StringIO()):
                s, M, Mf, P, cont, dec = A.assimilate_one(float(t), save_nc=False)
            if name == "v110":
                f = CONTROL_DIR / f"{t}Ma_corrected_SW.nc"
                with nc.Dataset(f) as d:
                    Mr = np.ma.filled(d["M_corrected"][:], np.nan).astype(float)
                dmax = float(np.nanmax(np.abs(Mr - Mf)))
                if not dmax <= 1.0:
                    raise SystemExit(f"[ea] CONTROL FAILED {t} Ma: v110 differs from {CONTROL_DIR.name} by {dmax:.2f} m")
            dla = np.asarray(lat[dec["iy"].astype(int).values]) if len(dec) else np.array([])
            dlo = np.asarray(lon[dec["ix"].astype(int).values]) if len(dec) else np.array([])

            def med3(la, lo):
                if not (np.isfinite(la) and np.isfinite(lo)):
                    return np.nan
                iy = int(nearest_cell_index(lat, np.array([la]))[0]); ix = int(nearest_cell_index(lon, np.array([lo]))[0])
                ys = slice(max(iy - 1, 0), iy + 2); xs = [(ix + d) % len(lon) for d in (-1, 0, 1)]
                return float(np.nanmedian(Mf[ys][:, xs]))

            # records at this age
            for i, r in recs.iterrows():
                if t not in r["_ages"]:
                    continue
                z = med3(site_lat[i], site_lon[i])
                dist = float(gc_km(site_lat[i], site_lon[i], dla, dlo).min()) if len(dla) and np.isfinite(site_lat[i]) else np.nan
                rec_rows.append(dict(config=name, age=t, rec=i, site=r["site"], use=r["use"], z=z, lo=r["lo_m"], hi=r["hi_m"],
                                     misfit=misfit(z, r["lo_m"], r["hi_m"]) if np.isfinite(z) else np.nan, dist_km=dist,
                                     prior=med3_prior(M, lat, lon, site_lat[i], site_lon[i], nearest_cell_index)))
            # East-West contrast and high bodies
            row = dict(config=name, age=t)
            for k in ("W", "E"):
                iy, ix = cells(*R[k]); sel = cont[iy, ix]
                row[f"median_{k}"] = float(np.median(Mf[iy[sel], ix[sel]])) if sel.any() else np.nan
            row["contrast_E_minus_W"] = row["median_E"] - row["median_W"]
            iy, ix = cells(*R["EA"])
            ea = np.zeros_like(cont, bool); ea[iy, ix] = True
            ea &= cont
            bb = bodies(ea & (Mf > HIGH_M))
            areas = sorted([float(cell_km2[b].sum()) for b in bb], reverse=True)
            row.update(high_bodies=len(bb), high_largest_km2=areas[0] if areas else 0.0, high_total_km2=float(sum(areas)),
                       ea_land_median=float(np.median(Mf[ea & (M > 0)])) if (ea & (M > 0)).any() else np.nan,
                       n_samples=int(len(dec)))
            age_rows.append(row)
            for j, bn in enumerate(BASINS):
                basin_rows.append(dict(config=name, age=t, basin=bn, z=med3(b_lat[j], b_lon[j])))
            print(f"[ea] {t:3d} Ma {name:17s} E-W {row['contrast_E_minus_W']:+6.0f} m; >3 km bodies {row['high_bodies']:2d}, "
                  f"largest {row['high_largest_km2'] / 1e3:6.0f} x10^3 km2", flush=True)
    return rec_rows, age_rows, basin_rows


def med3_prior(M, lat, lon, la, lo, nci):
    if not (np.isfinite(la) and np.isfinite(lo)):
        return np.nan
    iy = int(nci(lat, np.array([la]))[0]); ix = int(nci(lon, np.array([lo]))[0])
    ys = slice(max(iy - 1, 0), iy + 2); xs = [(ix + d) % len(lon) for d in (-1, 0, 1)]
    return float(np.nanmedian(M[ys][:, xs]))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--ages", default="200-65", help="'A-B' every grid age between, or a comma list")
    ap.add_argument("--configs", default=",".join(CONFIGS))
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--west", default=None, help="lat0,lat1,lon0,lon1 of the western box (present day); default WEST_BOX")
    ap.add_argument("--east", default=None, help="lat0,lat1,lon0,lon1 of the eastern box (present day); default EAST_BOX")
    ap.add_argument("--out", type=Path, default=OUTPUT_DIR / "east_asia_ablation")
    a = ap.parse_args()
    _self_test()
    if a.self_test:
        return 0
    global WEST_BOX, EAST_BOX
    if a.west:
        WEST_BOX = tuple(float(x) for x in a.west.split(","))
    if a.east:
        EAST_BOX = tuple(float(x) for x in a.east.split(","))
    print(f"[ea] boxes (present day, lat0, lat1, lon0, lon1): West {WEST_BOX}, East {EAST_BOX}")
    import pandas as pd
    sys.path.insert(0, str(HERE))
    from sw_io import available_ages
    grid = [int(round(t)) for t in available_ages() if t > 0]
    if "-" in a.ages:
        hi, lo = sorted(float(x) for x in a.ages.split("-"))[::-1]
        ages = [t for t in grid if lo <= t <= hi]
    else:
        ages = [int(x) for x in a.ages.split(",")]
    names = a.configs.split(",")
    if "v110" not in names:
        names = ["v110"] + names
    for n in names:
        if n not in CONFIGS:
            raise SystemExit(f"unknown configuration {n}")
    if a.out.exists() and any(a.out.iterdir()):
        raise SystemExit(f"[ea] {a.out} is not empty: move its contents elsewhere first")
    recs = pd.read_csv(RECORDS)
    for c in ("lo_m", "hi_m"):
        recs[c] = recs[c].astype(float)
    recs["_ages"] = [record_ages(r.age_min_ma, r.age_max_ma, ages) for r in recs.itertuples()]
    chunks = [ages[i::a.jobs] for i in range(a.jobs) if ages[i::a.jobs]]
    from concurrent.futures import ProcessPoolExecutor
    RR, AA, BB = [], [], []
    with ProcessPoolExecutor(len(chunks)) as ex:
        for r, s, b in ex.map(worker, [(c, names, recs, WEST_BOX, EAST_BOX) for c in chunks]):
            RR += r; AA += s; BB += b
    a.out.mkdir(parents=True, exist_ok=True)
    Rd, Ad, Bd = pd.DataFrame(RR), pd.DataFrame(AA).sort_values(["config", "age"]), pd.DataFrame(BB)
    for df, fn in ((Rd, "records.csv"), (Ad, "age_metrics.csv"), (Bd, "basins.csv")):
        tmp = a.out / (fn + ".partial"); df.to_csv(tmp, index=False, float_format="%.1f"); tmp.rename(a.out / fn)
    # R1: per record mean misfit over its ages, then mean over records, per configuration
    out = dict(west_box=WEST_BOX, east_box=EAST_BOX, rule="R1: v110 within 100 m of the better of rev_reach_support and submitted_fixed (strong East Asian records)",
               ages=ages, configs={})
    if len(Rd):
        per = Rd.groupby(["config", "rec", "use"]).agg(misfit=("misfit", "mean"), dist_km=("dist_km", "min")).reset_index()
        for name in names:
            p = per[per.config == name]
            st = p[p.use == "strong"]; far = st[st.dist_km > FAR_KM]
            out["configs"][name] = dict(strong_mean=float(st.misfit.mean()) if len(st) else None, strong_n=int(len(st)),
                                        strong_far_mean=float(far.misfit.mean()) if len(far) else None, strong_far_n=int(len(far)),
                                        weak_mean=float(p[p.use == "weak"].misfit.mean()) if (p.use == "weak").any() else None,
                                        tibet_mean=float(p[p.use == "tibet"].misfit.mean()) if (p.use == "tibet").any() else None,
                                        by_use={u: float(g.misfit.mean()) for u, g in p.groupby("use")})
        c = out["configs"]
        if all(k in c and c[k]["strong_mean"] is not None for k in ("v110", "rev_reach_support", "submitted_fixed")):
            best = min(c["rev_reach_support"]["strong_mean"], c["submitted_fixed"]["strong_mean"])
            out["R1"] = dict(v110=c["v110"]["strong_mean"], best_alternative=best,
                             verdict="KEEP the reach and support limits" if c["v110"]["strong_mean"] <= best + TOL_M
                             else "REVISE the far-field rule globally")
        for name in names:
            x = c.get(name, {})
            print(f"[ea] {name:17s} strong records: mean misfit {x.get('strong_mean')} m (n {x.get('strong_n')}); "
                  f"far (> {FAR_KM:.0f} km) {x.get('strong_far_mean')} (n {x.get('strong_far_n')}); weak {x.get('weak_mean')}; "
                  f"Tibet {x.get('tibet_mean')}")
        if "R1" in out:
            print(f"[ea] R1: v110 {out['R1']['v110']:.0f} m against the best alternative {out['R1']['best_alternative']:.0f} m "
                  f"-> {out['R1']['verdict']}")
    tmp = a.out / "summary.json.partial"; json.dump(out, open(tmp, "w"), indent=1, default=float); tmp.rename(a.out / "summary.json")
    print(f"[ea] wrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
