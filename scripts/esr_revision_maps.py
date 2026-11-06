#!/usr/bin/env python3
"""esr_revision_maps.py — maps of the input and corrected grids with their support layer (pyGMT, Winkel-Tripel).

  fig5      SW_comparison: six ages (500, 400, 300, 200, 100, 50 Ma), input, corrected and delta z, with
            - the support layer: on the corrected maps a white veil where the corrected grid is the input grid
              (> 500 km from the nearest declustered observation) and a lighter one where the residual kernel is
              tapered (250–500 km); on the delta-z maps the continental area is shaded in the same three classes
              beneath the correction;
            - continental outlines (the S&W18 continental mask, so flooded continent reads as continent) and the
              names of the principal blocks, placed by reconstructing a present-day point on each block's plate;
            - the declustered observations on the delta-z maps.

Inputs (paths from paths_scotese): the corrected grids in CORRECTED_DIR, whose support_km variable is the support
layer (OUTPUT_DIR/esr_revision/support_km.nc, written by esr_revision_diagnostics.py, supplies it for grids that do not
carry the variable, such as the v1.1.0 grids of the Zenodo archive), and OUTPUT_DIR/esr_revision/declustered.csv.
Output: FIGURES_DIR/SW_comparison.{png,pdf}, Fig. 5 of the paper (build_all_figures.sh gives it its number).

Label placement is checked for overlaps on the projected page before anything is drawn (a Winkel-Tripel forward
projection and an estimate of each label's box); the script stops if two labels would collide.

    cd scripts
    python3 esr_revision_maps.py fig5
    python3 esr_revision_maps.py --self-test
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from paths_scotese import CORRECTED_DIR as PRODUCT, OUTPUT_DIR, FIGURES_DIR as FIG  # noqa: E402

DIAG = OUTPUT_DIR / "esr_revision"
AGES = [500, 400, 300, 200, 100, 50]
ELEV_RANGE = (-4000.0, 4000.0, 250.0)
DELTA_RANGE = (-2000.0, 2000.0, 100.0)
SUPPORT_CLASSES = (250.0, 500.0)      # km: kernel whole / tapered / none (RESID_REACH_KM)
DZ_SHOWN_M = 20.0                     # |delta z| below this is not coloured on the delta-z maps

# Present-day anchor (lat, lon) on each block, its S&W plate ID (rotation-file comments), and the ages it is named at.
BLOCKS = [
    ("Laurentia", 45.0, -95.0, 101, (300, 600)), ("North America", 45.0, -100.0, 101, (0, 299)),
    ("Baltica", 62.0, 25.0, 301, (390, 600)), ("Siberia", 64.0, 105.0, 401, (150, 600)),
    ("Gondwana", 2.0, 20.0, 701, (300, 600)), ("Africa", 2.0, 20.0, 701, (0, 299)),
    ("South America", -12.0, -55.0, 201, (0, 299)), ("Eurasia", 52.0, 45.0, 301, (0, 389)),
    ("North China", 37.0, 112.0, 604, (150, 600)), ("South China", 27.0, 110.0, 611, (150, 600)),
    ("Kazakhstania", 47.0, 72.0, 402, (300, 600)), ("India", 20.0, 78.0, 501, (0, 200)),
    ("Australia", -25.0, 134.0, 801, (0, 200)), ("Antarctica", -80.0, 60.0, 802, (0, 200)),
]
LABEL_PT = 7.0


# ------------------------------------------------------------------------------------------------ geometry
def winkel(lon, lat, width_cm):
    """Winkel-Tripel (pygmt R0/<width>c) forward projection: page cm from the map's lower-left corner."""
    lon = (np.asarray(lon, float) + 180.0) % 360.0 - 180.0
    lam, phi = np.radians(lon), np.radians(np.asarray(lat, float))
    p1 = np.arccos(2.0 / np.pi)
    a = np.arccos(np.cos(phi) * np.cos(lam / 2.0))
    sa = np.sinc(a / np.pi)
    x = 0.5 * (lam * np.cos(p1) + 2.0 * np.cos(phi) * np.sin(lam / 2.0) / sa)
    y = 0.5 * (phi + np.sin(phi) / sa)
    xmax = 0.5 * np.pi * (np.cos(p1) + 1.0)                       # lon 180, lat 0
    ymax = np.pi / 2.0                                            # the pole
    s = width_cm / (2.0 * xmax)
    return s * (x + xmax), s * (y + ymax), 2.0 * xmax * s, 2.0 * ymax * s


def text_box(x, y, text, pt, justify):
    """Approximate page box (cm) of a Helvetica-Bold label: 0.60 em per character, 1.0 em tall."""
    w = 0.60 * pt * len(text) / 72.0 * 2.54
    h = 1.0 * pt / 72.0 * 2.54
    x0 = x - w / 2 if justify[1] == "C" else (x if justify[1] == "L" else x - w)
    return (x0, y - h / 2, x0 + w, y + h / 2)


def overlap(a, b, pad=0.05):
    return a[0] < b[2] + pad and b[0] < a[2] + pad and a[1] < b[3] + pad and b[1] < a[3] + pad


def place(labels, width_cm, fixed=()):
    """labels: list of (text, lon, lat). Returns [(text, lon, lat, dx_cm, dy_cm, justify)] with no two boxes
    overlapping each other or the fixed boxes; raises if a label cannot be placed."""
    out, boxes = [], list(fixed)
    offs = [(0, 0, "MC"), (0, 0.32, "MC"), (0, -0.32, "MC"), (0.15, 0, "ML"), (-0.15, 0, "MR"),
            (0, 0.6, "MC"), (0, -0.6, "MC"), (0.15, 0.32, "ML"), (-0.15, -0.32, "MR")]
    for text, lon, lat in labels:
        x, y, W, H = winkel(lon, lat, width_cm)
        for dx, dy, j in offs:
            b = text_box(x + dx, y + dy, text, LABEL_PT, j)
            if b[0] >= 0.05 and b[2] <= W - 0.05 and b[1] >= 0.05 and b[3] <= H - 0.05 and not any(overlap(b, o) for o in boxes):
                boxes.append(b); out.append((text, lon, lat, dx, dy, j)); break
        else:
            raise SystemExit(f"[esr-maps] no free place for the label {text!r} at ({lat:.0f}, {lon:.0f})")
    return out


def block_labels(t):
    import contextlib, io
    from sample_reconstruct_scotese import ScoteseSampleReconstructor
    with contextlib.redirect_stdout(io.StringIO()):
        sr = ScoteseSampleReconstructor()
    rows = [b for b in BLOCKS if b[4][0] <= t <= b[4][1]]
    df = pd.DataFrame(dict(Lat=[b[1] for b in rows], Lon=[b[2] for b in rows], sw_plate_id=[b[3] for b in rows]))
    rla, rlo = sr.reconstruct(df, float(t))
    return [(b[0], float(lo), float(la)) for b, la, lo in zip(rows, rla, rlo) if np.isfinite(la) and np.isfinite(lo)]


# ------------------------------------------------------------------------------------------------ data
def load(t):
    import netCDF4 as nc
    with nc.Dataset(PRODUCT / f"{t}Ma_corrected_SW.nc") as d:
        g = {k: np.ma.filled(d[k][:], np.nan).astype(float) for k in ("M_orig", "M_corrected", "delta")}
        g["cont"] = np.ma.filled(d["continent_mask"][:], 0) > 0
        g["lat"] = d["lat"][:].astype(float); g["lon"] = d["lon"][:].astype(float)
        if "support_km" in d.variables:
            g["support"] = np.ma.filled(d["support_km"][:], np.nan).astype(float)
    if "support" not in g:                     # grids without the support_km variable
        with nc.Dataset(DIAG / "support_km.nc") as d:
            k = int(np.where(np.round(d["age"][:]).astype(int) == t)[0][0])
            g["support"] = np.ma.filled(d["support_km"][k], np.nan).astype(float)
    dec = pd.read_csv(DIAG / "declustered.csv")
    g["dec"] = dec[dec.age == t]
    return g


def to_xr(arr, lat, lon):
    """As make_comparison_figures.to_xr: pixel-registered, geographic."""
    import xarray as xr
    da = xr.DataArray(arr.astype(np.float32), coords={"lat": lat, "lon": lon}, dims=("lat", "lon"), name="z")
    da.gmt.registration = 1
    da.gmt.gtype = 1
    return da


# ------------------------------------------------------------------------------------------------ Fig. 5
def fig5(width_cm=7.6, dpi=300):
    import pygmt
    FIG.mkdir(parents=True, exist_ok=True)
    tmp = FIG / "_tmp"; tmp.mkdir(exist_ok=True)
    elev_cpt, delta_cpt = str(tmp / "elev.cpt"), str(tmp / "delta.cpt")
    pygmt.makecpt(cmap="earth", series=ELEV_RANGE, continuous=True, background=True, output=elev_cpt)
    pygmt.makecpt(cmap="polar", series=DELTA_RANGE, continuous=True, background=True, output=delta_cpt)
    sup_cpt = str(tmp / "support.cpt")
    Path(sup_cpt).write_text("0 255/255/253 1 255/255/253\n1 224/224/221 2 224/224/221\n2 184/184/180 3 184/184/180\nN 127/127/127\n")
    # GMT 6.5 draws a CPT of pure greys as a one-channel image, and -Q (nan_transparent) then paints the NaN cells
    # red. A faint tint makes the image RGB and the NaN cells transparent.
    veil_cpt = str(tmp / "veil.cpt")
    Path(veil_cpt).write_text("0 255/255/253 3 255/255/253\nN 127/127/127\n")

    proj = f"R0/{width_cm}c"
    _, _, W, H = winkel(0, 0, width_cm)
    gap_x, gap_row, gap_group, cb_band = 0.45, 0.55, 2.3, 0.0
    fig = pygmt.Figure()
    pygmt.config(MAP_FRAME_TYPE="plain", MAP_FRAME_PEN="0.6p,black", FONT_ANNOT_PRIMARY="7p,Helvetica,black",
                 FONT_LABEL="8.5p,Helvetica,black", PS_PAGE_COLOR="white")
    cur = [0.0, 0.0]
    X0 = 0.6                                  # room for the row names left of the first column

    def goto(x, y):
        fig.shift_origin(xshift=f"{x - cur[0]:.3f}c", yshift=f"{y - cur[1]:.3f}c")
        cur[0], cur[1] = x, y

    letters = iter("abcdefghijklmnopqr")
    rows = ["Input", "Corrected", "@~D@~z"]
    stats = []
    # page layout: group 0 (500, 400, 300) on top, group 1 (200, 100, 50) below; within a group rows top to bottom
    Y0 = 2.6                                  # room for the colour bars and the key below the lower group
    y_top = Y0 + 2 * (3 * H + 2 * gap_row) + gap_group
    for gi, ages in enumerate((AGES[:3], AGES[3:])):
        for ci, t in enumerate(ages):
            g = load(t)
            lat, lon = g["lat"], g["lon"]
            cont = g["cont"]
            sup = g["support"]
            cls = np.where(sup <= SUPPORT_CLASSES[0], 0.5, np.where(sup <= SUPPORT_CLASSES[1], 1.5, 2.5))
            cls = np.where(np.isfinite(sup), cls, 2.5)
            cont_da = to_xr(cont.astype(np.float32), lat, lon)
            labels = place(block_labels(t), width_cm)
            dz = np.where(np.abs(g["delta"]) >= DZ_SHOWN_M, g["delta"], np.nan)
            w = np.cos(np.radians(lat))[:, None] * np.ones((1, len(lon)))
            stats.append((t, len(g["dec"]), *[float(((cont & (cls == k)) * w).sum() / (cont * w).sum())
                                              for k in (0.5, 2.5)]))
            for ri in range(3):
                x = X0 + ci * (W + gap_x)
                y = y_top - gi * (3 * H + 2 * gap_row + gap_group) - (ri + 1) * H - ri * gap_row
                goto(x, y)
                fig.basemap(region="g", projection=proj, frame=["f"])
                if ri == 0:
                    fig.grdimage(grid=to_xr(g["M_orig"], lat, lon), projection=proj, region="g", cmap=elev_cpt,
                                 nan_transparent=True)
                elif ri == 1:
                    fig.grdimage(grid=to_xr(g["M_corrected"], lat, lon), projection=proj, region="g", cmap=elev_cpt,
                                 nan_transparent=True)
                    for k, tr in ((1.5, 75), (2.5, 50)):              # tapered: light veil; beyond 500 km: strong
                        v = np.where(cont & (g["M_corrected"] >= 0) & (cls == k), 1.0, np.nan)   # on land only
                        fig.grdimage(grid=to_xr(v, lat, lon), projection=proj, region="g", cmap=veil_cpt,
                                     nan_transparent=True, transparency=tr)
                else:
                    fig.grdimage(grid=to_xr(np.where(cont, cls, np.nan), lat, lon), projection=proj, region="g",
                                 cmap=sup_cpt, nan_transparent=True)
                    fig.grdimage(grid=to_xr(dz, lat, lon), projection=proj, region="g", cmap=delta_cpt,
                                 nan_transparent=True)
                fig.grdcontour(grid=cont_da, projection=proj, region="g", levels=[0.5], pen="0.35p,gray20")
                if ri == 2 and len(g["dec"]):
                    fig.plot(x=g["dec"].cell_lon.values, y=g["dec"].cell_lat.values, style="c0.07c",
                             pen="0.4p,black", projection=proj, region="g")
                if ri == 0:
                    for text, lo, la, dx, dy, j in labels:
                        fig.text(x=lo, y=la, text=text, offset=f"{dx:.2f}c/{dy:.2f}c", justify=j,
                                 font=f"{LABEL_PT}p,Helvetica-Bold,black=~0.6p,white", region="g", projection=proj)
                fig.text(x=0.0, y=H + 0.08, text=next(letters), font="10p,Helvetica-Bold,black", justify="BL",
                         region=[0, W, 0, H], projection=f"X{W:.3f}c/{H:.3f}c", no_clip=True)
                if ri == 0:
                    fig.text(x=W / 2, y=H + 0.08, text=f"{t} Ma", font="10p,Helvetica-Bold,black", justify="BC",
                             region=[0, W, 0, H], projection=f"X{W:.3f}c/{H:.3f}c", no_clip=True)
                if ci == 0:
                    fig.text(x=-0.15, y=H / 2, text=rows[ri], font="9p,Helvetica-Bold,black", justify="BC",
                             angle=90, region=[0, W, 0, H], projection=f"X{W:.3f}c/{H:.3f}c", no_clip=True)
    # colour bars and the support key under the lower group
    goto(X0, Y0)
    tw = 3 * W + 2 * gap_x
    frame_box = dict(region=[0, tw, 0, 1], projection=f"X{tw:.3f}c/1c")
    fig.colorbar(cmap=elev_cpt, frame=["x+lelevation (m)", "af"], position=f"x0.2c/-0.9c+w{tw * 0.30:.2f}c/0.25c+h+e",
                 **frame_box)
    fig.colorbar(cmap=delta_cpt, frame=["x+l@~D@~z  corrected - input (m)", "af"],
                 position=f"x{tw * 0.36:.2f}c/-0.9c+w{tw * 0.30:.2f}c/0.25c+h+e", **frame_box)
    kx = tw * 0.72
    for i, (fill, lab) in enumerate((("255/255/253", "within 250 km: kernel whole"), ("224/224/221", "250-500 km: kernel tapered"),
                                     ("184/184/180", "beyond 500 km: input unchanged"))):
        fig.plot(x=[kx + 0.15], y=[-0.95 - 0.42 * i], style="s0.30c", fill=fill, pen="0.4p,black", no_clip=True,
                 **frame_box)
        fig.text(x=kx + 0.42, y=-0.95 - 0.42 * i, text=lab, font="7.5p,Helvetica,black", justify="ML", no_clip=True,
                 **frame_box)
    fig.text(x=kx, y=-0.55, text="distance to nearest declustered observation", font="7.5p,Helvetica-Bold,black",
             justify="ML", no_clip=True, **frame_box)
    out = FIG / "SW_comparison"
    fig.savefig(f"{out}.png", dpi=dpi)
    fig.savefig(f"{out}.pdf")
    for f in tmp.iterdir():
        f.unlink()
    tmp.rmdir()
    print(f"[esr-maps] wrote {out}.png/.pdf")
    for t, n, f0, f2 in stats:
        print(f"[esr-maps] {t:3d} Ma: {n:3d} declustered; continental area within 250 km {100 * f0:4.1f} %, "
              f"beyond 500 km {100 * f2:4.1f} %")


# ------------------------------------------------------------------------------------------------ self-test
def _self_test():
    from figure_checks import winkel_x_cm
    for lo, la in ((0, 0), (-180, 0), (120, 45), (-60, -70)):
        assert abs(winkel(lo, la, 10.0)[0] - winkel_x_cm(lo, la, 10.0)) < 1e-9
    x, y, W, H = winkel(0, 90, 10.0)
    assert abs(y - H) < 1e-9 and abs(x - 5.0) < 1e-9 and abs(W / H - 1.637) < 0.01, (W, H)
    p = place([("AAAA", 0, 0), ("BBBB", 0.5, 0)], 10.0)
    b0 = text_box(*winkel(0, 0, 10)[:2], "AAAA", LABEL_PT, p[0][5])
    x1, y1 = winkel(0.5, 0, 10)[:2]
    b1 = text_box(x1 + p[1][3], y1 + p[1][4], "BBBB", LABEL_PT, p[1][5])
    assert not overlap(b0, b1), "placed labels overlap"
    print("[esr-maps] self-test passed")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("what", nargs="*", default=["fig5"])
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    _self_test()
    if a.self_test:
        return 0
    for w in a.what:
        globals()[w]()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
