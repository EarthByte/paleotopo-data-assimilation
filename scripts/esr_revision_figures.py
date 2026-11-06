#!/usr/bin/env python3
"""esr_revision_figures.py — graphs and tables of the validation and of the corrected grids' support. The maps are in
esr_revision_maps.py, the regional figures in esr_revision_regions.py, the workflow diagram in esr_revision_flowchart.py.

Every number is read from files the workflow wrote (paths from paths_scotese):

  OUTPUT_DIR/esr_revision/          per_age.csv, declustered.csv (esr_revision_diagnostics.py);
                                    fig3_reconstructed_by_plate.csv (written by fig3 itself)
  OUTPUT_DIR/kernel_length_sweep/   reach_age_summary.csv, reach_heldout_samples.csv.gz: the held-out
                                    cross-validation (five folds of 6 degree blocks, every age 5-540 Ma) of the
                                    workflow with kernel lengths 150 km (candidate reach500_land, the product), 100,
                                    200 and 300 km (LS100, LS200, LS300); read by figS5 and tableS2 only
  data/independent_records/geology_record_comparison_50Myr_bins.csv   the Geology record (Zhou et al., 2026) on its
                                    50 Myr bins, with the area above 3 km of the corrected grids in the same bins
  CORRECTED_DIR                     the corrected grids

Figures (PNG 300 dpi + PDF) go to FIGURES_DIR under descriptive names; build_all_figures.sh gives them their paper
numbers. Tables (CSV + Markdown) go to OUTPUT_DIR/tables.

  fig3     Fig. 3, plate_id_validation: (a) samples by host plate; (b) per slice and host plate, the share of
           reconstructed samples inside the slice's continental mask; (c) share outside the mask and share on the input
           map's sea
  fig9     Fig. 14, area_above_thresholds: area above 2 and 3 km, input and corrected, log axis, lighter where < 30
           declustered observations; (b) declustered observations per slice
  fig10    Fig. 15, geology_record_comparison: the corrected grids against the Geology record on its 50 Myr bins
  figS5    Fig. S5 and Table S2, kernel_length_sweep (+ table): kernel length 100/150/200/300 km
  tableS2  Table S3 and Fig. S6, validation_by_province (+ figure): misfit by province, in sample and held out

    cd scripts
    python3 esr_revision_figures.py all          # or: fig3 fig9 fig10 figS5 tableS2
    python3 esr_revision_figures.py --self-test

No titles, panel letters only; every figure is checked for overlapping text before it is written.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from figure_checks import text_collisions, legend_over_data  # noqa: E402
from paths_scotese import PROJECT_ROOT as ROOT, CORRECTED_DIR as PRODUCT, OUTPUT_DIR, FIGURES_DIR as FIG  # noqa: E402

DIAG = OUTPUT_DIR / "esr_revision"
SWEEP = OUTPUT_DIR / "kernel_length_sweep"
GEOLOGY = ROOT / "data" / "independent_records" / "geology_record_comparison_50Myr_bins.csv"
TAB = OUTPUT_DIR / "tables"

ERAS = [(0, 66, "Cenozoic", "#fae6c8"), (66, 252, "Mesozoic", "#cfe6cf"), (252, 540, "Paleozoic", "#cfd9e6")]
WITH_SUBMITTED = False   # fig10: also draw the first version of the correction (pre_* columns; --with-submitted)
N_RELIABLE = 30          # declustered observations: the province rescaling's full strength (N_FULL_P)
CAND = {"LS100": 100, "reach500_land": 150, "LS200": 200, "LS300": 300}     # 150 km is the product
CAND_COL = {100: "#7fb3d5", 150: "#c0392b", 200: "#2e86c1", 300: "#1b4f72"}
PROV_ORDER = ["Continental Arc", "Island Arc", "Orogen", "Continental Margin", "Extended Crust", "Basin",
              "Platform", "Shield"]
plt.rcParams.update({"font.size": 10, "axes.labelsize": 10.5, "xtick.labelsize": 9.5, "ytick.labelsize": 9.5,
                     "legend.fontsize": 9, "pdf.fonttype": 42, "svg.fonttype": "none"})


# ------------------------------------------------------------------------------------------------ helpers
def letter(ax, s, x=-0.08, y=1.02):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=13, fontweight="bold", va="bottom", ha="right")


def eras(ax, names=True):
    for lo, hi, name, col in ERAS:
        ax.axvspan(lo, hi, color=col, alpha=0.55, zorder=0, lw=0)
        if names:
            ax.text((lo + hi) / 2, 1.005, name, transform=ax.get_xaxis_transform(), ha="center", va="bottom",
                    fontsize=9, color="#333")


def save(fig, name, legends=()):
    bad = text_collisions(fig)
    if bad:
        raise SystemExit(f"[esr-fig] {name}: overlapping text {bad[:6]}")
    for ax, leg in legends:
        n = legend_over_data(ax, leg)
        if n:
            raise SystemExit(f"[esr-fig] {name}: legend covers {n} plotted points")
    FIG.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"{name}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[esr-fig] wrote {FIG / name}.png/.pdf (text-collision check passed)")


def place_labels(ax, xs, ys, texts, avoid=None, **kw):
    """Label points without overlap: each label takes the first offset whose box clears every label placed before
    it and every labelled point; the figure check catches anything left."""
    fig = ax.figure
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    pts = ax.transData.transform(np.c_[xs, ys])
    if avoid is not None:
        pts = np.r_[pts, ax.transData.transform(np.asarray(avoid, float))]
    offs = [(5, -3, "left"), (5, 5, "left"), (-5, 5, "right"), (-5, -3, "right"), (5, -11, "left"),
            (-5, -11, "right"), (5, 13, "left"), (-5, 13, "right"), (12, -3, "left"), (-12, -3, "right")]
    placed = []
    for (x, y), t in zip(zip(xs, ys), texts):
        for dx, dy, ha in offs:
            a = ax.annotate(t, (x, y), xytext=(dx, dy), textcoords="offset points", ha=ha, va="center", **kw)
            bb = a.get_window_extent(r).expanded(1.05, 1.15)
            hit = any(bb.overlaps(p) for p in placed) or any(
                bb.x0 - 3 < px < bb.x1 + 3 and bb.y0 - 3 < py < bb.y1 + 3 for px, py in pts)
            if not hit:
                placed.append(bb)
                break
            a.remove()
        else:
            a = ax.annotate(t, (x, y), xytext=offs[0][:2], textcoords="offset points", va="center", **kw)
            placed.append(a.get_window_extent(r))


def rms(e):
    e = np.asarray(e, float)
    return float(np.sqrt(np.mean(e ** 2))) if len(e) else np.nan


def write_table(df, name, md_note=""):
    TAB.mkdir(parents=True, exist_ok=True)
    df.to_csv(TAB / f"{name}.csv", index=False, float_format="%.4g")
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "|".join("---" for _ in cols) + "|"]
    for _, r in df.iterrows():
        lines.append("| " + " | ".join(_fmt(r[c]) for c in cols) + " |")
    (TAB / f"{name}.md").write_text((md_note + "\n\n" if md_note else "") + "\n".join(lines) + "\n")
    print(f"[esr-fig] wrote {TAB / name}.csv/.md")


def _fmt(v):
    if isinstance(v, (float, np.floating)):
        if not np.isfinite(v):
            return "–"
        return f"{v:,.0f}" if abs(v) >= 100 else (f"{v:.2f}" if abs(v) >= 0.01 or v == 0 else f"{v:.3f}")
    return str(v)


# ------------------------------------------------------------------------------------------------ Fig. 9
def fig9():
    d = pd.read_csv(DIAG / "per_age.csv").sort_values("age")
    floor = 0.003                                   # % ; zero area is drawn on the floor with open markers
    fig, (a, b) = plt.subplots(2, 1, figsize=(8.6, 6.4), sharex=True, gridspec_kw=dict(height_ratios=[2.3, 1]))
    eras(a)
    ok = d.n_decluster >= N_RELIABLE
    series = [("frac_gt2_corrected", "#555555", "-", 1.8, "corrected, > 2 km"),
              ("frac_gt3_corrected", "#c0392b", "-", 2.0, "corrected, > 3 km"),
              ("frac_gt2_input", "#555555", "--", 1.0, "input, > 2 km"),
              ("frac_gt3_input", "#c0392b", "--", 1.0, "input, > 3 km")]
    handles = []
    for col, c, ls, lw, lab in series:
        y = 100 * d[col].values
        yp = np.where(y > 0, y, np.nan)
        a.plot(d.age, yp, color=c, ls=ls, lw=lw, alpha=0.35, zorder=2)                 # whole series, light
        h, = a.plot(d.age, np.where(ok, yp, np.nan), color=c, ls=ls, lw=lw, zorder=3, label=lab)   # n >= 30, full
        handles.append(h)
        z = (y == 0) & col.endswith("corrected")
        if z.any():
            a.plot(d.age[z], np.full(z.sum(), floor * 1.25), ls="none", marker="o", mfc="none", mec=c, ms=3.5,
                   zorder=3)
    a.set_yscale("log")
    a.set_ylim(floor, 12)
    a.set_ylabel("Continental area (%)")
    a.grid(True, which="major", alpha=0.3, lw=0.6)
    leg = a.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False, borderaxespad=0)
    a.text(1.02, 0.30, f"lighter: fewer than {N_RELIABLE}\ndeclustered observations\n\nopen circles: no corrected\narea above 3 km",
           transform=a.transAxes, fontsize=8.5, va="top", color="#333")
    letter(a, "a")
    eras(b, names=False)
    b.bar(d.age, d.n_decluster, width=4.2, color=np.where(ok, "#34495e", "#aab7c4"), zorder=2)
    b.axhline(N_RELIABLE, color="#c0392b", lw=1.0, ls="--", zorder=3)
    b.text(1.006, N_RELIABLE, f"{N_RELIABLE}", color="#c0392b", fontsize=9, va="center", ha="left",
           transform=b.get_yaxis_transform())
    b.set_yscale("log")
    b.set_ylim(8, 1000)
    b.set_ylabel("Declustered\nobservations")
    b.set_xlabel("Age (Ma)")
    b.set_xlim(540, 0)
    b.grid(True, which="major", axis="y", alpha=0.3, lw=0.6)
    letter(b, "b")
    fig.subplots_adjust(hspace=0.08)
    save(fig, "area_above_thresholds")


# ------------------------------------------------------------------------------------------------ Fig. 10
def geology_stats(g):
    from scipy.stats import pearsonr, spearmanr
    out = {}
    for k in ("v110", "pre"):
        r, p = pearsonr(g.geo_p95_km, g[f"{k}_frac_gt3"])
        rho, prho = spearmanr(g.geo_p95_km, g[f"{k}_frac_gt3"])
        top_geo = set(g.nlargest(4, "geo_p95_km").bin_lo)
        top_k = set(g.nlargest(4, f"{k}_frac_gt3").bin_lo)
        out[k] = dict(r=r, p=p, rho=rho, prho=prho, shared_top4=len(top_geo & top_k), top_k=sorted(top_k))
    return out


def fig10(with_submitted=None):
    """The paper's Fig. 15: the released grids against the sample record. --with-submitted also draws the first version
    of the correction (the pre_* columns of the record file)."""
    if with_submitted is None:
        with_submitted = WITH_SUBMITTED
    g = pd.read_csv(GEOLOGY).sort_values("bin_lo")
    st = geology_stats(g)
    mid = (g.bin_lo + g.bin_hi) / 2
    fig = plt.figure(figsize=(10.2, 5.6))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.55, 1], hspace=0.10, wspace=0.30)
    a = fig.add_subplot(gs[0, 0]); b = fig.add_subplot(gs[1, 0], sharex=a); c = fig.add_subplot(gs[:, 1])
    for ax in (a, b):
        eras(ax, names=ax is a)
        for _, r in g[g.superhighland].iterrows():
            ax.axvspan(r.bin_lo, r.bin_hi, ymin=0, ymax=1, facecolor="none", edgecolor="#7d3c98", hatch="////",
                       lw=0, alpha=0.35, zorder=1)
    edges = np.r_[g.bin_lo.values, g.bin_hi.values[-1]]
    a.stairs(g.geo_p95_km, edges, color="#7d3c98", lw=2.0, zorder=3, baseline=None, label="95th percentile")
    a.stairs(g.geo_median_km, edges, color="#7d3c98", lw=1.2, ls="--", zorder=3, baseline=None, label="median")
    from matplotlib.patches import Patch
    a.set_ylabel("Sample elevation (km)")
    a.set_ylim(0, 6.2)
    hh, ll = a.get_legend_handles_labels()
    hh.append(Patch(facecolor="none", edgecolor="#b39bc8", hatch="////", lw=0)); ll.append("super-highland bin")
    leg_a = a.legend(hh, ll, loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=3, frameon=True, framealpha=0.92,
                     edgecolor="none", fontsize=8.5, columnspacing=1.2, handlelength=1.8)
    a.tick_params(labelbottom=False)
    a.grid(True, axis="y", alpha=0.3, lw=0.6)
    letter(a, "a")
    w = 18
    if with_submitted:
        b.bar(mid - w / 2, g.v110_frac_gt3, width=w, color="#c0392b", zorder=3, label="corrected")
        b.bar(mid + w / 2, g.pre_frac_gt3, width=w, color="#e6b0aa", zorder=3, label="submitted")
    else:
        b.bar(mid, g.v110_frac_gt3, width=2 * w, color="#c0392b", zorder=3, label="corrected grids")
    b.set_yscale("log")
    b.set_ylim(0.01, 5)
    b.set_ylabel("Area > 3 km (%)")
    b.set_xlabel("Age (Ma)")
    b.set_xlim(550, 0)
    b.grid(True, axis="y", which="major", alpha=0.3, lw=0.6)
    leg_b = b.legend(loc="upper left", frameon=False, ncol=2)
    letter(b, "b")
    # (c) scatter
    series = (("pre", "o", "#e6b0aa", "submitted"), ("v110", "o", "#c0392b", "corrected"))
    for k, mk, col, lab in (series if with_submitted else series[1:]):
        c.scatter(g.geo_p95_km, g[f"{k}_frac_gt3"], s=46 if k == "v110" else 30, marker=mk,
                  facecolor=col if k == "v110" else "none", edgecolor=col, lw=1.3, zorder=3, label=lab)
    c.set_yscale("log")
    c.set_xlim(3.85, 5.45)
    c.set_ylim(0.015, 3.0)
    place_labels(c, g.geo_p95_km.values, g.v110_frac_gt3.values,
                 [f"{int(r.bin_lo)}–{int(r.bin_hi)}" for _, r in g.iterrows()],
                 avoid=np.c_[g.geo_p95_km, g.pre_frac_gt3] if with_submitted else None,
                 fontsize=7.5, color="#444")
    c.set_yscale("log")
    c.set_xlim(3.85, 5.45)
    c.set_ylim(0.015, 3.0)
    c.set_xlabel("Sample elevation, 95th percentile (km)")
    c.set_ylabel("Corrected area > 3 km (%)")
    s = st["v110"]
    txt = f"r = {s['r']:+.2f} (p = {s['p']:.3f})\nρ = {s['rho']:+.2f} (p = {s['prho']:.2f})"
    if with_submitted:
        txt = f"corrected: {txt}\nsubmitted: r = {st['pre']['r']:+.2f}, ρ = {st['pre']['rho']:+.2f}"
    c.text(0.03, 0.97, txt, transform=c.transAxes, va="top", fontsize=9.5, family="DejaVu Sans")
    leg_c = c.legend(loc="lower right", frameon=False) if with_submitted else None
    c.grid(True, which="major", alpha=0.3, lw=0.6)
    letter(c, "c", x=-0.14)
    save(fig, "geology_record_comparison" + ("_with_submitted" if with_submitted else ""),
         legends=[(c, leg_c)] if leg_c is not None else [])
    return st


# ------------------------------------------------------------------------------------------------ Fig. S5 / Table S3
def figS5():
    s = pd.read_csv(SWEEP / "reach_age_summary.csv")
    h = pd.read_csv(SWEEP / "reach_heldout_samples.csv.gz")
    s["ls"] = s.candidate.map(CAND); h["ls"] = h.candidate.map(CAND)
    s, h = s.dropna(subset=["ls"]), h.dropna(subset=["ls"])
    L = h[h.land_sited]
    bins = [0, 150, 250, 400, 600, 20000]
    blab = ["0–150", "150–250", "250–400", "400–600", "> 600"]
    L = L.assign(db=pd.cut(L.dist_land_km, bins, right=True, labels=blab, include_lowest=True))
    fig = plt.figure(figsize=(10.4, 7.2))
    gs = fig.add_gridspec(2, 2, hspace=0.32, wspace=0.28)
    a = fig.add_subplot(gs[0, 0]); b = fig.add_subplot(gs[0, 1]); c = fig.add_subplot(gs[1, :])
    x = np.arange(len(blab))
    ref = L[L.ls == 150]
    pr = [rms(ref[ref.db == lb].z - ref[ref.db == lb].prior) for lb in blab]
    a.plot(x, pr, color="#888", lw=1.6, ls="--", marker="s", ms=4, label="input (prior)")
    for ls in sorted(CAND.values()):
        q = L[L.ls == ls]
        a.plot(x, [rms(q[q.db == lb].z - q[q.db == lb].pred) for lb in blab], color=CAND_COL[ls], lw=1.8,
               marker="o", ms=4, label=f"{ls} km" + (" (product)" if ls == 150 else ""))
    nn = ref.groupby("db", observed=False).size()
    a.set_xticks(x, [f"{lb}\nn = {nn[lb]:,}" for lb in blab], fontsize=8.5)
    a.set_xlabel("Distance to nearest training sample (km)")
    a.set_ylabel("Held-out RMS misfit (m)")
    a.set_ylim(1000, 2300)
    a.grid(True, alpha=0.3, lw=0.6)
    leg_a = a.legend(loc="lower right", frameon=False, fontsize=8.5)
    letter(a, "a", x=-0.13)
    # (b) in-sample vs held-out, medians over ages
    for ls in sorted(CAND.values()):
        q = s[s.ls == ls]
        b.scatter(q.insample_rms.median(), np.median(q.heldout_land_near_rms.dropna()), s=70, color=CAND_COL[ls],
                  zorder=3, label=f"{ls} km")
        b.annotate(f"{ls} km", (q.insample_rms.median(), np.median(q.heldout_land_near_rms.dropna())),
                   xytext=(7, 3), textcoords="offset points", fontsize=8.5)
    b.set_xlabel("In-sample RMS, median over ages (m)")
    b.set_ylabel("Held-out RMS ≤ 150 km,\nmedian over ages (m)")
    b.grid(True, alpha=0.3, lw=0.6)
    xl = b.get_xlim(); b.set_xlim(xl[0] - 0.15 * (xl[1] - xl[0]), xl[1] + 0.35 * (xl[1] - xl[0]))
    yl = b.get_ylim(); b.set_ylim(yl[0] - 0.25 * (yl[1] - yl[0]), yl[1] + 0.25 * (yl[1] - yl[0]))
    letter(b, "b", x=-0.16)
    # (c) area above 3 km through time per length scale
    eras(c)
    for ls in sorted(CAND.values()):
        q = s[s.ls == ls].sort_values("age")
        y = 100 * q.frac_gt3_after.values
        c.plot(q.age, np.where(y > 0, y, np.nan), color=CAND_COL[ls], lw=2.0 if ls == 150 else 1.2,
               label=f"{ls} km", zorder=3 if ls == 150 else 2)
    c.set_yscale("log")
    c.set_ylim(0.003, 12)
    c.set_xlim(540, 0)
    c.set_xlabel("Age (Ma)")
    c.set_ylabel("Land area > 3 km (%)")
    c.grid(True, which="major", alpha=0.3, lw=0.6)
    leg_c = c.legend(loc="upper left", bbox_to_anchor=(1.005, 1.0), frameon=False, borderaxespad=0)
    letter(c, "c", x=-0.06)
    save(fig, "kernel_length_sweep")
    # Table S3
    rows = []
    for ls in sorted(CAND.values()):
        q, Lq = s[s.ls == ls], L[L.ls == ls]
        rows.append({"kernel length (km)": ls,
                     "held-out n (land)": len(Lq),
                     "held-out RMS, all (m)": rms(Lq.z - Lq.pred),
                     "held-out RMS ≤ 150 km (m)": rms(Lq[Lq.dist_land_km <= 150].z - Lq[Lq.dist_land_km <= 150].pred),
                     "held-out RMS > 400 km (m)": rms(Lq[Lq.dist_land_km > 400].z - Lq[Lq.dist_land_km > 400].pred),
                     "in-sample RMS, median (m)": q.insample_rms.median(),
                     "land p99, median (m)": q.land_p99_after.median(),
                     "area > 2 km, median (%)": 100 * q.frac_gt2_after.median(),
                     "area > 3 km, median (%)": 100 * q.frac_gt3_after.median(),
                     "lowland lifted > 300 m, median (%)": 100 * q.lowland_lifted_300.median()})
    L150 = L[L.ls == 150]
    rows.append({"kernel length (km)": "input (prior)", "held-out n (land)": len(L150),
                 "held-out RMS, all (m)": rms(L150.z - L150.prior),
                 "held-out RMS ≤ 150 km (m)": rms(L150[L150.dist_land_km <= 150].z - L150[L150.dist_land_km <= 150].prior),
                 "held-out RMS > 400 km (m)": rms(L150[L150.dist_land_km > 400].z - L150[L150.dist_land_km > 400].prior),
                 "in-sample RMS, median (m)": s[s.ls == 150].prior_rms.median(),
                 "land p99, median (m)": s[s.ls == 150].land_p99_before.median(),
                 "area > 2 km, median (%)": np.nan, "area > 3 km, median (%)": np.nan,
                 "lowland lifted > 300 m, median (%)": 0.0})
    write_table(pd.DataFrame(rows), "kernel_length_sweep_table",
                "Kernel length sweep, every age 5–540 Ma (108 ages), the product's workflow otherwise unchanged. "
                "Held-out: five folds of 6° blocks, samples on the input map's land; distance to the nearest training "
                "sample on land. Areas are fractions of land area.")


# ------------------------------------------------------------------------------------------------ Table S2 / Fig. S6
def tableS2():
    h = pd.read_csv(SWEEP / "reach_heldout_samples.csv.gz")
    h = h[(h.candidate == "reach500_land") & h.land_sited]
    dc = pd.read_csv(DIAG / "declustered.csv")
    rows = []
    for p in PROV_ORDER + ["All"]:
        I = dc if p == "All" else dc[dc.prov == p]
        H = h if p == "All" else h[h.prov == p]
        N = H[H.dist_land_km <= 150]
        if len(I) != len(H):
            raise SystemExit(f"[esr-fig] tableS2 {p}: {len(I)} declustered against {len(H)} held-out land samples; "
                             "the two sets should be the same observations")
        rows.append({"province": p, "n": len(I),
                     "input bias (m)": float((I.prior - I.z).mean()),
                     "input RMS (m)": rms(I.prior - I.z),
                     "in-sample bias (m)": float((I.corrected - I.z).mean()),
                     "in-sample RMS (m)": rms(I.corrected - I.z),
                     "held-out bias (m)": float((H.pred - H.z).mean()),
                     "held-out RMS (m)": rms(H.pred - H.z),
                     "≤ 150 km n": len(N),
                     "≤ 150 km input RMS (m)": rms(N.prior - N.z),
                     "≤ 150 km held-out RMS (m)": rms(N.pred - N.z)})
    T = pd.DataFrame(rows)
    write_table(T, "validation_by_province_table",
                "Validation by tectonic province (150 km kernel). Bias = model − observation. "
                "The same 7,894 declustered observations (all slices 5–540 Ma) are scored three ways: against the "
                "input grid; in sample, against the corrected grid built with them; held out, against the corrected "
                "grid built without them (five folds of 6° blocks). ≤ 150 km: held-out "
                "observations within 150 km of the nearest training observation on land.")
    # Fig. S6: the same as a dot chart
    T = T[T.province != "All"]
    fig, axs = plt.subplots(1, 3, figsize=(10.6, 3.9), sharey=True)
    y = np.arange(len(T))[::-1]
    for ax, (k0, k1, kn, lab), let in zip(axs, (("input RMS (m)", "in-sample RMS (m)", "n", "In sample"),
                                                ("input RMS (m)", "held-out RMS (m)", "n", "Held out, all"),
                                                ("≤ 150 km input RMS (m)", "≤ 150 km held-out RMS (m)", "≤ 150 km n",
                                                 "Held out, ≤ 150 km")), "abc"):
        i0, i1 = T[k0].values, T[k1].values
        for yy, u, v in zip(y, i0, i1):
            ax.annotate("", xy=(v, yy), xytext=(u, yy), arrowprops=dict(arrowstyle="-|>", color="#c0392b", lw=1.3))
        ax.scatter(i0, y, s=28, facecolor="white", edgecolor="#555", zorder=3, label="input")
        ax.scatter(i1, y, s=30, color="#c0392b", zorder=3, label="corrected")
        n = T[kn].values
        for yy, v, nn in zip(y, np.nanmax(np.c_[i0, i1], 1), n):
            ax.text(v + 60, yy, f"{nn:,}", va="center", fontsize=8, color="#555")
        ax.set_xlabel(f"{lab}: RMS misfit (m)")
        ax.grid(True, axis="x", alpha=0.3, lw=0.6)
        ax.set_xlim(0, 2500); ax.set_xticks([0, 500, 1000, 1500, 2000])
        letter(ax, let, x=-0.03 if let != "a" else -0.05)
    axs[0].set_yticks(y, T.province)
    axs[2].legend(loc="lower left", bbox_to_anchor=(1.0, 0.0), frameon=False, fontsize=8.5, borderaxespad=0.3)
    fig.subplots_adjust(wspace=0.10)
    save(fig, "validation_by_province")


# ------------------------------------------------------------------------------------------------ Fig. 3
# Names from the comments of Scotese_Wright_PlateModel_2023.rot (604 NCH North China, 611 SCH South China, 613 LHA Lhasa,
# 628 AMU Amuria; 615 is Indochina).
PLATE_NAMES = {101: "North America", 201: "South America", 291: "Paraná", 301: "Eurasia", 302: "Baltica",
               304: "Iberia", 305: "Variscan Europe", 401: "Siberia", 402: "Kazakhstania", 501: "India",
               601: "Tarim", 602: "Qaidam", 604: "North China", 611: "South China", 612: "Qiangtang", 613: "Lhasa",
               615: "Indochina", 628: "Amuria", 701: "Africa", 801: "Australia", 802: "East Antarctica"}


def fig3_data():
    """Per slice and host plate: reconstructed samples, inside the slice's continental mask, on the input map's land."""
    import contextlib, io
    import netCDF4 as nc
    import assimilate_scotese as A
    from sw_io import nearest_cell_index, available_ages
    with contextlib.redirect_stdout(io.StringIO()):
        df = A.get_geochem()
    rows = []
    for t in [int(round(x)) for x in available_ages()]:
        if t == 0:
            continue
        with contextlib.redirect_stdout(io.StringIO()):
            s = A.prepare_samples(df, t)
        if s.empty:
            continue
        with nc.Dataset(PRODUCT / f"{t}Ma_corrected_SW.nc") as d:
            cont = np.ma.filled(d["continent_mask"][:], 0) > 0
            M = np.ma.filled(d["M_orig"][:], np.nan).astype(float)
            lat, lon = d["lat"][:], d["lon"][:]
        iy = nearest_cell_index(lat, s.rlat.values); ix = nearest_cell_index(lon, s.rlon.values)
        g = pd.DataFrame(dict(pid=s.sw_plate_id.astype(int).values, inside=cont[iy, ix], land=M[iy, ix] > 0))
        for pid, q in g.groupby("pid"):
            rows.append(dict(age=t, pid=pid, n=len(q), inside=int(q.inside.sum()), land=int(q.land.sum())))
    out = pd.DataFrame(rows)
    DIAG.mkdir(parents=True, exist_ok=True)
    out.to_csv(DIAG / "fig3_reconstructed_by_plate.csv", index=False)
    return out, df


def fig3():
    import cartopy.crs as ccrs
    import matplotlib.colors as mcolors
    f = DIAG / "fig3_reconstructed_by_plate.csv"
    import contextlib, io
    d, df = fig3_data()
    tot = d.groupby("pid").n.sum().sort_values(ascending=False)
    top = tot.index[:12].tolist()
    fig = plt.figure(figsize=(10.6, 10.4))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.45, 1.15, 0.80], hspace=0.16)
    # (a) present-day samples by plate
    a = fig.add_subplot(gs[0], projection=ccrs.Robinson())
    a.set_global(); a.coastlines(lw=0.4, color="grey")
    pc = df.sw_plate_id.value_counts(); major = pc.head(8).index.tolist()
    pal = [plt.cm.tab10(k) for k in (0, 1, 2, 3, 4, 5, 6, 8, 9)]   # tab10 without its grey
    for i, pid in enumerate(major + [-1]):
        sub = df[df.sw_plate_id.isin(major) == False] if pid == -1 else df[df.sw_plate_id == pid]  # noqa: E712
        col = (0.72, 0.72, 0.72, 0.8) if pid == -1 else pal[i]
        a.scatter(sub.Lon, sub.Lat, s=3, color=col, alpha=0.6, edgecolor="none", transform=ccrs.PlateCarree(),
                  label="other plates" if pid == -1 else f"{int(pid)} {PLATE_NAMES.get(int(pid), '')}".strip())
    leg_a = a.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False, markerscale=4, fontsize=8.5,
                     borderaxespad=0)
    for hnd in leg_a.legend_handles:
        hnd.set_alpha(1.0)
    letter(a, "a", x=-0.02)
    # (b) heat map: share inside the continental mask, per host plate and slice
    b = fig.add_subplot(gs[1])
    ages = sorted(d.age.unique())
    Z = np.full((len(top), len(ages)), np.nan)
    for i, pid in enumerate(top):
        q = d[d.pid == pid].set_index("age")
        for j, t in enumerate(ages):
            if t in q.index and q.loc[t, "n"] >= 3:
                Z[i, j] = 100 * q.loc[t, "inside"] / q.loc[t, "n"]
    edges = np.r_[np.array(ages) - 2.5, ages[-1] + 2.5]
    cm = mcolors.LinearSegmentedColormap.from_list("ins", ["#7b241c", "#e59866", "#fdebd0", "#d5f5e3", "#1e8449"])
    im = b.pcolormesh(edges, np.arange(len(top) + 1), Z, cmap=cm, vmin=50, vmax=100, shading="flat")
    nu = df.sw_plate_id.value_counts()
    b.set_yticks(np.arange(len(top)) + 0.5, [f"{p} {PLATE_NAMES.get(p, '')}".strip() + f"  ({nu.get(p, 0):,})" for p in top],
                 fontsize=8.5)
    b.invert_yaxis()
    b.set_xlim(540, 0)
    b.set_xlabel("Age (Ma)")
    cax = b.inset_axes([1.015, 0.0, 0.018, 1.0])
    cb = fig.colorbar(im, cax=cax, extend="min"); cb.set_label("Inside continental mask (%)", fontsize=9)
    b.text(1.0, -0.035, "white: fewer than 3 samples in the slice", transform=b.transAxes, ha="right", va="top",
           fontsize=8.5, color="#444")
    letter(b, "b", x=-0.30)
    # (c) all samples: outside the mask, on the input map's sea
    c = fig.add_subplot(gs[2], sharex=b)
    eras(c, names=False)
    A_ = d.groupby("age")[["n", "inside", "land"]].sum()
    c.plot(A_.index, 100 * (1 - A_.inside / A_.n), color="#7b241c", lw=1.6, label="outside the continental mask")
    c.plot(A_.index, 100 * (1 - A_.land / A_.n), color="#2e86c1", lw=1.6, label="on the input map's sea")
    c.axhline(25, color="grey", lw=0.8, ls="--", label="25 % (submitted warning level)")
    c.set_ylim(0, 100)
    c.set_ylabel("Samples (%)")
    c.set_xlabel("Age (Ma)")
    c.set_xlim(540, 0)
    c.grid(True, alpha=0.3, lw=0.6)
    leg_c = c.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False, borderaxespad=0, fontsize=8.5)
    letter(c, "c", x=-0.30)
    b.tick_params(labelbottom=False); b.set_xlabel("")
    save(fig, "plate_id_validation")
    s = A_.assign(out=100 * (1 - A_.inside / A_.n), sea=100 * (1 - A_.land / A_.n))
    print(f"[esr-fig] fig3: outside mask median {s.out.median():.1f} % (max {s.out.max():.1f} % at {s.out.idxmax()} Ma); "
          f"on sea median {s.sea.median():.1f} %")


# ------------------------------------------------------------------------------------------------ self-test
def _self_test():
    assert abs(rms([3, 4]) - np.sqrt(12.5)) < 1e-12 and np.isnan(rms([]))
    g = pd.DataFrame(dict(bin_lo=[0, 50, 100, 150, 200], bin_hi=[50, 100, 150, 200, 250],
                          geo_p95_km=[5, 4, 3, 2, 1], v110_frac_gt3=[1, .5, .2, .1, .05], pre_frac_gt3=[1, .1, .2, .3, .4]))
    s = geology_stats(g)
    assert abs(s["v110"]["rho"] - 1) < 1e-12 and s["v110"]["shared_top4"] == 4, s
    fig, ax = plt.subplots(); ax.text(.1, .5, "AAAA"); ax.text(.11, .5, "BBBB")
    assert text_collisions(fig); plt.close(fig)
    print("[esr-fig] self-test passed")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("what", nargs="*", default=["all"])
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--with-submitted", action="store_true",
                    help="fig10: also draw the grids as first submitted (the response to the reviewers)")
    a = ap.parse_args()
    global WITH_SUBMITTED
    WITH_SUBMITTED = a.with_submitted
    _self_test()
    if a.self_test:
        return 0
    todo = ["fig9", "fig10", "figS5", "tableS2", "fig3"] if "all" in a.what else a.what
    for w in todo:
        r = globals()[w]()
        if w == "fig10":
            for k, v in r.items():
                print(f"[esr-fig] fig10 {k}: r {v['r']:+.2f} (p {v['p']:.3f}), rho {v['rho']:+.2f} (p {v['prho']:.2f}); "
                      f"top-4 bins {v['top_k']}, {v['shared_top4']} shared with the record's top 4")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
