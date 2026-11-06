#!/usr/bin/env python3
"""esr_revision_regions.py — regional figures: well-known orogenic belts and arcs, input against corrected, with the
support layer, the samples and independent paleoelevation records.

One figure per region. Rows: three ages. Columns: input (S&W18), corrected (with the support veil and the declustered
observations coloured by their own elevation, on the same colour scale), Δz (on the three support classes).
Independent paleoelevation records of the region are numbered on the corrected map of the age nearest their midpoint,
and the bottom panel compares every record with the input and corrected grids at its reconstructed position (3 x 3
cell median, mean over the grid ages inside the record's age interval; as east_asia_ablation.py).

  andes        central Andes since 60 Ma
  tibet        southern Tibet since 60 Ma
  cordillera   North American Cordillera, Late Cretaceous
  cpm          Central Pangean Mountains (Variscan-Appalachian), 320-280 Ma
  east_asia    eastern Asia, Jurassic-Cretaceous

Records: data/independent_records/<region>.csv (schema in data/independent_records/README.md). Only methods independent
of whole-rock igneous geochemistry are admitted, because igneous geochemistry is what the correction is trained on.
A region without its file is drawn without records.

Output (paths from paths_scotese): FIGURES_DIR/region_<name>.{png,pdf}, and the record-by-record comparison in
OUTPUT_DIR/tables/region_<name>_records.csv. In the paper these are Figs. 8-12 and Tables S4-S8, in the order andes,
tibet, cordillera, cpm, east_asia.

    cd scripts
    python3 esr_revision_regions.py all        # or: andes tibet cordillera cpm east_asia
"""
from __future__ import annotations

import argparse
import contextlib
import io
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from esr_revision_maps import load, to_xr, SUPPORT_CLASSES, DZ_SHOWN_M, ELEV_RANGE, DELTA_RANGE   # noqa: E402
from paths_scotese import PROJECT_ROOT as ROOT, OUTPUT_DIR, FIGURES_DIR as FIG                   # noqa: E402

TAB = OUTPUT_DIR / "tables"
REC = ROOT / "data" / "independent_records"

# centre: a present-day point on the region's plate (S&W plate ID from the rotation-file comments), reconstructed to
# each age; radius: the map's horizon in degrees.
REGIONS = {
    "andes": dict(ages=[50, 30, 10], centre=(-18.0, -67.0, 201), radius=16, records="central_andes_paleoelevation_60-0Ma.csv",
                  zlim=(-500, 6000)),
    "tibet": dict(ages=[55, 35, 15], centre=(31.0, 88.0, 613), radius=17, records="tibet_paleoelevation_60-0Ma.csv",
                  zlim=(-500, 8000)),
    "cordillera": dict(ages=[95, 80, 70], centre=(40.0, -114.0, 101), radius=18,
                       records="cordillera_paleoelevation_100-66Ma.csv", zlim=(-500, 5000)),
    "cpm": dict(ages=[320, 300, 280], centre=(47.0, -10.0, 304), radius=30,
                records="central_pangean_mountains_paleoelevation_330-270Ma.csv", zlim=(-500, 5000)),
    "east_asia": dict(ages=[160, 130, 100], centre=(40.0, 118.0, 604), radius=20,
                      records="east_asia_paleoelevation_200-66Ma.csv", zlim=(-500, 6000)),
}


def misfit(z, lo, hi):
    return float(max(lo - z, 0.0, z - hi))


def record_ages(amin, amax, grid_ages):
    """As east_asia_ablation.record_ages: grid ages inside [amin, amax], else the nearest within 5 Myr of the midpoint."""
    a = [t for t in grid_ages if amin - 1e-9 <= t <= amax + 1e-9]
    if a:
        return a
    mid = 0.5 * (amin + amax)
    near = min(grid_ages, key=lambda t: abs(t - mid)) if grid_ages else None
    return [near] if near is not None and abs(near - mid) <= 5.0 else []


def med3(F, lat, lon, la, lo):
    from sw_io import nearest_cell_index
    if not (np.isfinite(la) and np.isfinite(lo)):
        return np.nan
    iy = int(nearest_cell_index(lat, np.array([la]))[0]); ix = int(nearest_cell_index(lon, np.array([lo]))[0])
    ys = slice(max(iy - 1, 0), iy + 2); xs = [(ix + d) % len(lon) for d in (-1, 0, 1)]
    return float(np.nanmedian(F[ys][:, xs]))


def aeqd_cm(lat, lon, clat, clon, radius_deg, width_cm):
    """Azimuthal equidistant page position (cm from the map's lower-left corner) as GMT -JE<clon>/<clat>/<radius>/<width>."""
    p1, l1, p2, l2 = map(np.radians, (clat, clon, np.asarray(lat, float), np.asarray(lon, float)))
    c = np.arccos(np.clip(np.sin(p1) * np.sin(p2) + np.cos(p1) * np.cos(p2) * np.cos(l2 - l1), -1, 1))
    az = np.arctan2(np.sin(l2 - l1) * np.cos(p2), np.cos(p1) * np.sin(p2) - np.sin(p1) * np.cos(p2) * np.cos(l2 - l1))
    r = width_cm / 2.0 * np.degrees(c) / radius_deg
    return width_cm / 2.0 + r * np.sin(az), width_cm / 2.0 + r * np.cos(az)


def record_labels(here, clat, clon, radius, W, pt=6.0, join_cm=0.45):
    """Cluster record sites closer than join_cm on the page into one label ("7,8,10"), then give each label the first
    offset whose box clears the other labels. Returns rows (x_cm, y_cm, label_x_cm, label_y_cm, text). Raises if a
    label finds no free place."""
    x, y = aeqd_cm(here.la.values, here.lo.values, clat, clon, radius, W)
    inside = np.hypot(x - W / 2, y - W / 2) <= 0.47 * W          # beyond the horizon: scored in the record panel only
    cl = []
    for xi, yi, n in sorted(zip(x[inside], y[inside], here.n.values[inside]), key=lambda q: q[2]):
        for c in cl:
            if np.hypot(xi - c[0], yi - c[1]) < join_cm:
                c[2].append(int(n)); break
        else:
            cl.append([xi, yi, [int(n)]])
    boxes, out = [], []
    offs = [(0, 0.32), (0, -0.32), (0.42, 0), (-0.42, 0), (0.36, 0.30), (-0.36, 0.30), (0.36, -0.30), (-0.36, -0.30),
            (0, 0.62), (0, -0.62)]
    for cx, cy, ns in cl:
        txt = ",".join(str(k) for k in ns)
        w = 0.60 * pt * len(txt) / 72 * 2.54 + 0.12; h = pt / 72 * 2.54 + 0.10
        for dx, dy in offs:
            b = (cx + dx - w / 2, cy + dy - h / 2, cx + dx + w / 2, cy + dy + h / 2)
            if not any(b[0] < o[2] + 0.03 and o[0] < b[2] + 0.03 and b[1] < o[3] + 0.03 and o[1] < b[3] + 0.03 for o in boxes):
                boxes.append(b); out.append((cx, cy, cx + dx, cy + dy, txt)); break
        else:
            raise SystemExit(f"[esr-regions] no free place for the record label {txt}")
    return out


def reconstructor():
    import assimilate_scotese as A
    with contextlib.redirect_stdout(io.StringIO()):
        return A.get_reconstructor()


def score_records(recs, sr):
    """Per record: input and corrected grid at the reconstructed site, mean over the grid ages of its interval."""
    from sw_io import available_ages
    grid_ages = [int(round(t)) for t in available_ages()]
    df = pd.DataFrame({"Lat": recs.lat.values, "Lon": recs.lon.values})
    df["sw_plate_id"] = sr.assign_plate_ids(df)
    out = []
    cache = {}
    for i, r in recs.reset_index(drop=True).iterrows():
        ages = record_ages(r.age_min_ma, r.age_max_ma, grid_ages)
        zi, zc, pos = [], [], []
        for t in ages:
            if t not in cache:
                cache[t] = load(t)
            g = cache[t]
            la, lo = sr.reconstruct(df.iloc[[i]], float(t))
            zi.append(med3(g["M_orig"], g["lat"], g["lon"], la[0], lo[0]))
            zc.append(med3(g["M_corrected"], g["lat"], g["lon"], la[0], lo[0]))
            pos.append((la[0], lo[0]))
        zi_m = float(np.nanmean(zi)) if len(zi) and np.isfinite(zi).any() else np.nan
        zc_m = float(np.nanmean(zc)) if len(zc) and np.isfinite(zc).any() else np.nan
        out.append(dict(n=i + 1, site=r.site, age_min_ma=r.age_min_ma, age_max_ma=r.age_max_ma, lo_m=r.lo_m, hi_m=r.hi_m,
                        use=r.get("use", ""), method=r.method, reference=r.reference_short, grid_ages=len(ages),
                        input_m=zi_m, corrected_m=zc_m,
                        misfit_input_m=misfit(zi_m, r.lo_m, r.hi_m) if np.isfinite(zi_m) else np.nan,
                        misfit_corrected_m=misfit(zc_m, r.lo_m, r.hi_m) if np.isfinite(zc_m) else np.nan))
    return pd.DataFrame(out)


def region(name, width_cm=6.2, dpi=300):
    import pygmt
    cfg = REGIONS[name]
    sr = reconstructor()
    f = REC / cfg["records"]
    recs = pd.read_csv(f) if f.exists() else None
    if recs is not None and name == "east_asia":
        recs = recs[recs.use != "tibet"]               # the Tibetan rows of that file are Cretaceous, outside Fig. Tibet
    if recs is not None:
        # rows without a numeric interval ("NEAR-MODERN:" or "RELATIVE:" in the notes) cannot be scored and are left
        # out; a single age stands for the interval
        recs = recs[np.isfinite(recs.age_min_ma) & (np.isfinite(recs.lo_m) | np.isfinite(recs.hi_m))].copy()
        recs["age_max_ma"] = recs.age_max_ma.fillna(recs.age_min_ma)
        recs = recs.reset_index(drop=True)
    sc = score_records(recs, sr) if recs is not None and len(recs) else None
    if sc is not None:
        TAB.mkdir(parents=True, exist_ok=True)
        sc.to_csv(TAB / f"region_{name}_records.csv", index=False, float_format="%.0f")
    FIG.mkdir(parents=True, exist_ok=True)
    tmp = FIG / f"_tmp_{name}"; tmp.mkdir(exist_ok=True)
    elev_cpt, delta_cpt = str(tmp / "elev.cpt"), str(tmp / "delta.cpt")
    pygmt.makecpt(cmap="earth", series=ELEV_RANGE, continuous=True, background=True, output=elev_cpt)
    pygmt.makecpt(cmap="polar", series=DELTA_RANGE, continuous=True, background=True, output=delta_cpt)
    sup_cpt, veil_cpt = str(tmp / "support.cpt"), str(tmp / "veil.cpt")
    # tinted, not pure grey: see esr_revision_maps.fig5 (GMT 6.5 -Q on a one-channel image)
    Path(sup_cpt).write_text("0 255/255/253 1 255/255/253\n1 224/224/221 2 224/224/221\n2 184/184/180 3 184/184/180\nN 127/127/127\n")
    Path(veil_cpt).write_text("0 255/255/253 3 255/255/253\nN 127/127/127\n")
    W = width_cm; gx, gy = 0.35, 0.45
    fig = pygmt.Figure()
    pygmt.config(MAP_FRAME_TYPE="plain", MAP_FRAME_PEN="0.6p,black", FONT_ANNOT_PRIMARY="7.5p,Helvetica,black",
                 FONT_LABEL="8.5p,Helvetica,black", MAP_GRID_PEN_PRIMARY="0.25p,gray45")
    X0, Y0 = 0.9, 8.6 if sc is not None else 2.4
    cur = [0.0, 0.0]

    def goto(x, y):
        fig.shift_origin(xshift=f"{x - cur[0]:.3f}c", yshift=f"{y - cur[1]:.3f}c")
        cur[0], cur[1] = x, y
    letters = iter("abcdefghijklmnop")
    rec_df = pd.DataFrame({"Lat": recs.lat.values, "Lon": recs.lon.values}) if sc is not None else None
    if rec_df is not None:
        rec_df["sw_plate_id"] = sr.assign_plate_ids(rec_df)
    for ri, t in enumerate(cfg["ages"]):
        g = load(t)
        lat, lon, cont = g["lat"], g["lon"], g["cont"]
        c = pd.DataFrame({"Lat": [cfg["centre"][0]], "Lon": [cfg["centre"][1]], "sw_plate_id": [cfg["centre"][2]]})
        cla, clo = sr.reconstruct(c, float(t))
        proj = f"E{clo[0]:.2f}/{cla[0]:.2f}/{cfg['radius']}/{W}c"
        sup = g["support"]
        cls = np.where(sup <= SUPPORT_CLASSES[0], 0.5, np.where(sup <= SUPPORT_CLASSES[1], 1.5, 2.5))
        cls = np.where(np.isfinite(sup), cls, 2.5)
        dz = np.where(np.abs(g["delta"]) >= DZ_SHOWN_M, g["delta"], np.nan)
        cont_da = to_xr(cont.astype(np.float32), lat, lon)
        y = Y0 + (len(cfg["ages"]) - 1 - ri) * (W + gy)
        here = None
        if sc is not None:                               # records whose interval is nearest this row's age
            mids = 0.5 * (recs.age_min_ma.values + recs.age_max_ma.values)
            row_of = np.array([int(np.argmin([abs(m - a) for a in cfg["ages"]])) for m in mids])
            idx = np.where(row_of == ri)[0]
            if len(idx):
                la, lo = sr.reconstruct(rec_df.iloc[idx], float(t))
                here = pd.DataFrame(dict(n=idx + 1, la=la, lo=lo)).dropna()
        for ci in range(3):
            goto(X0 + ci * (W + gx), y)
            fig.basemap(region="g", projection=proj, frame=["f"])
            if ci == 0:
                fig.grdimage(grid=to_xr(g["M_orig"], lat, lon), projection=proj, region="g", cmap=elev_cpt)
            elif ci == 1:
                fig.grdimage(grid=to_xr(g["M_corrected"], lat, lon), projection=proj, region="g", cmap=elev_cpt)
                for k, tr in ((1.5, 75), (2.5, 50)):
                    v = np.where(cont & (g["M_corrected"] >= 0) & (cls == k), 1.0, np.nan)
                    fig.grdimage(grid=to_xr(v, lat, lon), projection=proj, region="g", cmap=veil_cpt,
                                 nan_transparent=True, transparency=tr)
            else:
                fig.grdimage(grid=to_xr(np.where(cont, cls, np.nan), lat, lon), projection=proj, region="g",
                             cmap=sup_cpt, nan_transparent=True)
                fig.grdimage(grid=to_xr(dz, lat, lon), projection=proj, region="g", cmap=delta_cpt, nan_transparent=True)
            fig.grdcontour(grid=cont_da, projection=proj, region="g", levels=[0.5], pen="0.4p,gray20")
            fig.basemap(region="g", projection=proj, frame=["f"])   # no graticule: GMT 6.5 draws g10 on -JE with region g as dense stripes
            if ci == 1 and len(g["dec"]):
                d = g["dec"]
                fig.plot(x=d.cell_lon.values, y=d.cell_lat.values, fill=d.z.values, cmap=elev_cpt, style="c0.16c",
                         pen="0.4p,black", projection=proj, region="g")
            if ci == 2 and len(g["dec"]):
                fig.plot(x=g["dec"].cell_lon.values, y=g["dec"].cell_lat.values, style="c0.09c", pen="0.4p,black",
                         projection=proj, region="g")
            box = dict(region=[0, W, 0, W], projection=f"X{W:.3f}c/{W:.3f}c", no_clip=True)
            if ci == 1 and here is not None and len(here):
                for sx, sy, lx, ly, txt in record_labels(here, cla[0], clo[0], cfg["radius"], W):
                    fig.plot(x=[sx, lx], y=[sy, ly], pen="0.5p,black", **box)
                    fig.plot(x=[sx], y=[sy], style="s0.14c", fill="white", pen="0.6p,black", **box)
                    fig.text(x=lx, y=ly, text=txt, font="6p,Helvetica-Bold,black", fill="white", pen="0.4p,black",
                             clearance="0.04c/0.03c", **box)
            fig.text(x=0.0, y=W + 0.08, text=next(letters), font="10p,Helvetica-Bold,black", justify="BL", **box)
            if ri == 0:
                fig.text(x=W / 2, y=W + 0.08, text=["Input (S&W18)", "Corrected", "@~D@~z"][ci],
                         font="9p,Helvetica-Bold,black", justify="BC", **box)
            if ci == 0:
                fig.text(x=-0.2, y=W / 2, text=f"{t} Ma", font="10p,Helvetica-Bold,black", justify="BC", angle=90, **box)
    # colour bars
    tw = 3 * W + 2 * gx
    goto(X0, Y0)
    fb = dict(region=[0, tw, 0, 1], projection=f"X{tw:.3f}c/1c")
    fig.colorbar(cmap=elev_cpt, frame=["x+lelevation (m); circles: declustered observations", "af"],
                 position=f"x0.0c/-0.95c+w{tw * 0.42:.2f}c/0.22c+h+e", **fb)
    fig.colorbar(cmap=delta_cpt, frame=["x+l@~D@~z (m); grey: farther than 250 / 500 km", "af"],
                 position=f"x{tw * 0.55:.2f}c/-0.95c+w{tw * 0.42:.2f}c/0.22c+h+e", **fb)
    # record panel
    if sc is not None:
        s = sc.dropna(subset=["corrected_m"]).copy()
        ph, pw = 5.2, tw - 1.4
        goto(X0 + 1.2, 1.2)
        z0, z1 = cfg["zlim"]
        reg = [0.4, len(sc) + 0.6, z0, z1]
        fig.basemap(region=reg, projection=f"X{pw:.2f}c/{ph}c",
                    frame=["WSne", "xa1+lindependent record (numbered on the maps)", "yaf+lelevation (m)"])
        for _, r in sc.iterrows():
            fig.plot(x=[r.n, r.n], y=[r.lo_m, min(r.hi_m, z1)], pen="5p,gray80", region=reg, projection=f"X{pw:.2f}c/{ph}c")
        fig.plot(x=s.n.values, y=s.input_m.values, style="c0.20c", pen="0.8p,gray30", region=reg, projection=f"X{pw:.2f}c/{ph}c")
        fig.plot(x=s.n.values, y=s.corrected_m.values, style="c0.20c", fill="red3", pen="0.4p,black", region=reg,
                 projection=f"X{pw:.2f}c/{ph}c")
        # the key (bars, open and red circles) and the mean misfits are in the caption: the figures carry no titles,
        # and a header line here would sit on the top-axis ticks
        fig.text(x=reg[0], y=z1, text=next(letters), font="10p,Helvetica-Bold,black", justify="BR", offset="-0.9c/0.1c",
                 region=reg, projection=f"X{pw:.2f}c/{ph}c", no_clip=True)
    out = FIG / f"region_{name}"
    fig.savefig(f"{out}.png", dpi=dpi)
    fig.savefig(f"{out}.pdf")
    for q in tmp.iterdir():
        q.unlink()
    tmp.rmdir()
    print(f"[esr-regions] wrote {out}.png/.pdf" + (f"; {len(sc)} records, mean misfit input {np.nanmean(sc.misfit_input_m):.0f} m, "
                                                 f"corrected {np.nanmean(sc.misfit_corrected_m):.0f} m" if sc is not None else
                                                 "; no record file"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("what", nargs="*", default=["all"])
    a = ap.parse_args()
    for w in (list(REGIONS) if "all" in a.what else a.what):
        region(w)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
