#!/usr/bin/env python3
"""esr_revision_flowchart.py — the workflow diagram, as vector graphics (PDF, SVG, PNG).

Inputs, one-time preprocessing, the per-slice loop of eight steps and the outputs. The S&W18 model has resolved plate
boundaries for 0-100 Ma only, used solely for the continental-arc province, and the input box says so. The step texts
and the footer read their parameters from assimilate_scotese.py, and the script refuses to draw if those are not the
settings of the product of record. Layout is computed: every box is sized to its text, connectors are orthogonal, and
the figure is checked for overlapping boxes, overlapping text and text spilling out of its box before it is written.

    cd scripts
    python3 esr_revision_flowchart.py          # -> FIGURES_DIR/methodology_flowchart.{pdf,svg,png} (Fig. 2)
"""
from __future__ import annotations

import contextlib
import io
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from paths_scotese import FIGURES_DIR as OUT  # noqa: E402
with contextlib.redirect_stdout(io.StringIO()):
    import assimilate_scotese as A
from figure_checks import text_collisions

assert A.DROP_CRATON_SAMPLES and A.DROP_SEA_SITED_SAMPLES and A.DELTA_SMOOTH_SIGMA_FINAL == 0 and \
    tuple(A.RESID_REACH_KM) == (250.0, 500.0) and A.RESID_LS_KM == 150.0 and A.LOCAL_SUPPORT_KM == 400.0, \
    "assimilate_scotese.py does not hold the settings of the product of record"

plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42, "svg.fonttype": "none"})
C = {"input": "#dcdcdc", "stage": "#cfd9e6", "iter": "#ffe4b5", "output": "#cfe6cf"}
FS = 8.6                       # box text, pt
LINE = 0.168                   # inch per text line at FS, with leading
PAD = 0.16                     # inch, vertical padding inside a box
GAP = 0.17                     # inch between stacked boxes
W_FIG = 8.5                    # inch
COLS = {"in": (0.10, 2.55), "mid": (2.95, 5.55), "out": (5.95, 8.40)}   # x0, x1 in inch

INPUTS = [
    "Geochemical sample database\n(~100k igneous rocks; major and\ntrace elements, present-day position,\nage, tectonic province)",
    "Plate model: rotation file and\npresent-day static continental polygons;\nresolved plate boundaries for 0–100 Ma\nonly, used for the continental-arc\nprovince",
    "Kinematic paleotopography\n(Scotese and Wright 2018 PaleoDEMs,\nor another series with its\nplate model)",
]
PRE = [
    "Crustal-thickness mohometry\n(Zhou et al. 2025, machine learning\non 66 geochemical features)",
    "Airy isostasy, six-model ensemble:\nper-sample elevation z ± σ",
    "Plate-ID assignment (point in polygon,\nstatic continental polygons,\nanchor plate 000)",
    "Drop samples outside every polygon\n(~2.5 %) and inside craton outlines\n(Shirmard et al. 2025)",
]
LOOP = "For each target age t = 0, 5, 10, …, 540 Ma"
STEPS = [
    "1. Sample window |Age − t| ≤ Δt(t),\ntriangular weight in |Age − t|",
    "2. Reconstruct sample positions to t\n(plate ID and rotation model)",
    "3. Drop samples outside the slice's\ncontinental mask",
    "4. Decluster (1° × province bins,\nweighted median, Kish count); drop\nbins on the input map's sea",
    f"5. Tectonic-province grid\n(1-cell speckle removed)",
    f"6. Province-wise CDF rescaling\n(N_eff {A.N_MIN_P}→{A.N_FULL_P}), zero beyond {A.LOCAL_SUPPORT_KM:.0f} km\nof the nearest sample; province-edge\n"
    f"taper, cross-boundary σ {A.CROSS_BOUNDARY_SIGMA:g} cells",
    f"7. Residual kernel ({A.RESID_LS_KM:.0f} km Gaussian),\nwhole within {A.RESID_REACH_KM[0]:.0f} km of the nearest\nsample, zero at {A.RESID_REACH_KM[1]:.0f} km",
    f"8. Continent-masked Δz smoothing\n(σ {A.DELTA_SMOOTH_SIGMA:g} cells), cap |Δz| ≤ {A.DZ_MAX / 1000:g} km;\ncorrection ramped to zero at 0 Ma\nover the youngest {A.YOUNG_TAPER_T_FULL:g} Myr",
]
OUTPUTS = [
    "Corrected NetCDF per slice\n(input and corrected elevation,\nΔz, provinces, continental mask,\nn_eff, support_km)",
    "Per-slice statistics\n(bias and RMS before and after,\nhypsometric percentiles)",
    "Validation: held-out samples,\nkernel-length sweep,\nindependent paleoelevation\nrecords",
    "Figures and movies\n(input, corrected, Δz, samples,\nsupport layer)",
    "Applications: area above\nelevation thresholds, hypsometry,\nderived crustal thickness",
]


def wrap(text, width, heading=True):
    """Re-wrap a box text to the column: a newline after a box's first line (its heading) is kept, the rest flows."""
    import textwrap
    for a in ("σ ", "N_eff ", "et al. ", "Zhou et ", "1-cell ", "|Δz| ", "≤ ", "Shirmard et "):
        text = text.replace(a, a[:-1] + "\u00a0")
    if not heading:
        return "\n".join(textwrap.wrap(text.replace("\n", " "), width))
    head, _, rest = text.partition("\n")
    out = textwrap.wrap(head, width)
    if rest:
        out += textwrap.wrap(rest.replace("\n", " "), width)
    return "\n".join(out)


def h_of(text):
    return PAD * 2 + LINE * (text.count("\n") + 1)


def main():
    global INPUTS, PRE, STEPS, OUTPUTS
    INPUTS = [wrap(t, 33) for t in INPUTS]; OUTPUTS = [wrap(t, 33) for t in OUTPUTS]
    PRE = [wrap(t, 38, False) for t in PRE]; STEPS = [wrap(t, 38, False) for t in STEPS]
    ys_in = [h_of(t) for t in INPUTS]; ys_pre = [h_of(t) for t in PRE]; ys_st = [h_of(t) for t in STEPS]
    ys_out = [h_of(t) for t in OUTPUTS]
    loop_h = 0.42
    H_mid = sum(ys_pre) + GAP * len(ys_pre) + loop_h + GAP + sum(ys_st) + GAP * (len(ys_st) - 1)
    top_pad, bot_pad = 0.45, 0.55
    H_FIG = H_mid + top_pad + bot_pad
    fig = plt.figure(figsize=(W_FIG, H_FIG))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W_FIG); ax.set_ylim(0, H_FIG); ax.set_axis_off()
    boxes = []

    def box(x0, x1, ytop, h, text, col, weight="normal"):
        p = FancyBboxPatch((x0, ytop - h), x1 - x0, h, boxstyle="round,pad=0,rounding_size=0.06",
                           facecolor=col, edgecolor="black", lw=0.8)
        ax.add_patch(p)
        ax.text((x0 + x1) / 2, ytop - h / 2, text, ha="center", va="center", fontsize=FS, weight=weight,
                linespacing=1.25)
        boxes.append((x0, ytop - h, x1, ytop, text.split("\n")[0]))
        return (x0, ytop - h, x1, ytop)

    def arrow(xy0, xy1, col="#444", lw=1.0):
        ax.annotate("", xy=xy1, xytext=xy0, arrowprops=dict(arrowstyle="-|>", color=col, lw=lw,
                                                            shrinkA=0, shrinkB=0, mutation_scale=9))

    def line(xs, ys, col="#444", lw=1.0):
        ax.plot(xs, ys, color=col, lw=lw, solid_capstyle="butt")

    y_top = H_FIG - top_pad
    for name, x in (("INPUTS", COLS["in"]), ("PREPROCESSING AND ASSIMILATION", COLS["mid"]), ("OUTPUTS", COLS["out"])):
        ax.text((x[0] + x[1]) / 2, y_top + 0.16, name, ha="center", va="bottom", fontsize=8.6, weight="bold",
                color="#555")
    # middle column
    mx0, mx1 = COLS["mid"]; mxc = (mx0 + mx1) / 2
    y = y_top; mid = []
    for t, h in zip(PRE, ys_pre):
        mid.append(box(mx0, mx1, y, h, t, C["stage"])); y -= h + GAP
    loop = box(mx0 - 0.25, mx1 + 0.25, y, loop_h, LOOP, C["iter"], weight="bold"); y -= loop_h + GAP
    steps = []
    for t, h in zip(STEPS, ys_st):
        steps.append(box(mx0, mx1, y, h, t, C["stage"])); y -= h + GAP
    chain = mid + [loop] + steps
    for a, b in zip(chain[:-1], chain[1:]):
        arrow((mxc, a[1]), (mxc, b[3]))
    # inputs, aligned with the box they feed
    ix0, ix1 = COLS["in"]
    targets = [mid[0], mid[2], steps[0]]
    ins = []
    yi = y_top
    for t, h, tg in zip(INPUTS, ys_in, targets):
        ytop = min(yi, (tg[1] + tg[3]) / 2 + h / 2)
        ins.append(box(ix0, ix1, ytop, h, t, C["input"])); yi = ytop - h - GAP
    xm = (ix1 + mx0) / 2
    for k, (b, tg) in enumerate(zip(ins, targets)):
        yb, yt = (b[1] + b[3]) / 2, (tg[1] + tg[3]) / 2
        x_el = ix1 + 0.10 + 0.06 * k
        line([ix1, x_el, x_el], [yb, yb, yt]); arrow((x_el, yt), (mx0, yt))
    # plate model also feeds step 2 (reconstruction): a branch down the same elbow
    yb2 = (steps[1][1] + steps[1][3]) / 2
    x_el = ix1 + 0.16
    line([x_el, x_el], [(mid[2][1] + mid[2][3]) / 2, yb2]); arrow((x_el, yb2), (mx0, yb2))
    # loop-back on the left of the steps
    xl = mx0 - 0.10
    line([steps[-1][0], xl, xl], [(steps[-1][1] + steps[-1][3]) / 2] * 2 + [loop[1] - 0.0], col="#aa6e10")
    ax.annotate("", xy=(xl, loop[1]), xytext=(xl, loop[1] - 0.25),
                arrowprops=dict(arrowstyle="-|>", color="#aa6e10", lw=1.0, shrinkA=0, shrinkB=0, mutation_scale=9))
    ax.text(xl - 0.07, (steps[0][3] + steps[-1][1]) / 2, "next slice", rotation=90, ha="right", va="center",
            fontsize=8.2, color="#aa6e10")
    # outputs: stacked beside the steps, fed by one bus from step 8
    ox0, ox1 = COLS["out"]
    span_top, span_bot = steps[0][3], steps[-1][1]
    tot = sum(ys_out); g = (span_top - span_bot - tot) / (len(ys_out) - 1)
    yo = span_top; outs = []
    for t, h in zip(OUTPUTS, ys_out):
        outs.append(box(ox0, ox1, yo, h, t, C["output"])); yo -= h + g
    xb = (mx1 + ox0) / 2
    y8 = (steps[-1][1] + steps[-1][3]) / 2
    line([mx1, xb], [y8, y8], col="#1a6f1a")
    yc_all = [(o[1] + o[3]) / 2 for o in outs] + [y8]
    line([xb, xb], [min(yc_all), max(yc_all)], col="#1a6f1a")
    for o in outs:
        yc = (o[1] + o[3]) / 2
        arrow((xb, yc), (ox0, yc), col="#1a6f1a")
    # footer: parameters of the product
    ax.text(W_FIG / 2, 0.26,
            f"Product parameters: Δt(t) = 5/10/20 Myr;  N_min = {A.N_MIN_P}, N_full = {A.N_FULL_P};\n"
            f"rescaling support {A.LOCAL_SUPPORT_KM:.0f} km;  residual kernel {A.RESID_LS_KM:.0f} km, reach "
            f"{A.RESID_REACH_KM[0]:.0f}–{A.RESID_REACH_KM[1]:.0f} km;  |Δz| ≤ {A.DZ_MAX / 1000:g} km",
            ha="center", va="center", fontsize=7.8, style="italic", color="#333")
    # checks: boxes must not overlap one another, text must not overlap
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            a, b = boxes[i], boxes[j]
            if a[0] < b[2] - 1e-6 and b[0] < a[2] - 1e-6 and a[1] < b[3] - 1e-6 and b[1] < a[3] - 1e-6:
                raise SystemExit(f"[esr-flow] boxes overlap: {a[4]!r} and {b[4]!r}")
    bad = text_collisions(fig)
    if bad:
        raise SystemExit(f"[esr-flow] overlapping text {bad[:6]}")
    r = fig.canvas.get_renderer()
    for t in ax.texts:                    # every box text must fit inside its box
        bb = t.get_window_extent(r).transformed(ax.transData.inverted())
        for b in boxes:
            if b[0] <= (bb.x0 + bb.x1) / 2 <= b[2] and b[1] <= (bb.y0 + bb.y1) / 2 <= b[3]:
                if bb.x0 < b[0] + 0.03 or bb.x1 > b[2] - 0.03 or bb.y0 < b[1] + 0.02 or bb.y1 > b[3] - 0.02:
                    fig.savefig(OUT / "_debug_flowchart.png", dpi=110)
                    raise SystemExit(f"[esr-flow] text spills out of its box: {b[4]!r}")
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "svg"):
        fig.savefig(OUT / f"methodology_flowchart.{ext}")
    fig.savefig(OUT / "methodology_flowchart.png", dpi=300)
    print(f"[esr-flow] wrote {OUT}/methodology_flowchart.pdf/.svg/.png (checks passed)")


if __name__ == "__main__":
    main()
