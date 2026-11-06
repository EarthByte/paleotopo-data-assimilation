#!/usr/bin/env python3
"""figure_checks.py — checks the figure scripts run on their own output before they report success.

    text_collisions(fig)          matplotlib: pairs of overlapping text (tick labels, axis labels, legend, annotations)
    legend_over_data(ax, leg)     matplotlib: number of plotted points that fall under the legend box
    edge_clipped(png)             dark pixels on the outer rim of a saved PNG (for figures whose artwork never reaches
                                  the edge; NOT for pygmt maps, whose frame and ocean fill legitimately touch the crop)
    winkel_x_cm(lon, lat, width)  x position (cm from the map's left edge) of a point on a Winkel-Tripel map centred on
                                  0 degrees, the projection pygmt draws with R0/<width>c
    label_left_cm(...)            where the left end of a panel label lands on the page; below 0 it is cut off, because
                                  GMT drops what lies left of the page origin

    python3 scripts/figure_checks.py --self-test
"""
from __future__ import annotations

import argparse
import sys

import numpy as np


def text_collisions(fig):
    """Pairs of overlapping visible text boxes on a matplotlib figure."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    boxes = []
    for ax in fig.axes:
        # tick labels whose tick lies outside the axis view are never drawn, though matplotlib still lists them
        x0, x1 = sorted(ax.get_xlim()); y0, y1 = sorted(ax.get_ylim())
        xt = [t for t in ax.get_xticklabels() if x0 - 1e-9 <= t.get_position()[0] <= x1 + 1e-9]
        yt = [t for t in ax.get_yticklabels() if y0 - 1e-9 <= t.get_position()[1] <= y1 + 1e-9]
        items = xt + yt + [ax.xaxis.label, ax.yaxis.label] + list(ax.texts)
        leg = ax.get_legend()
        if leg is not None:
            items += list(leg.get_texts()) + [leg.get_title()]
        for t in items:
            if t.get_visible() and t.get_text().strip():
                boxes.append((t.get_text(), t.get_window_extent(r)))
    bad = []
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            a, b = boxes[i][1], boxes[j][1]
            if a.x0 < b.x1 - 0.5 and b.x0 < a.x1 - 0.5 and a.y0 < b.y1 - 0.5 and b.y0 < a.y1 - 0.5:
                bad.append((boxes[i][0], boxes[j][0]))
    return bad


def legend_over_data(ax, leg):
    """Number of plotted line vertices (display coordinates) inside the legend's box."""
    ax.figure.canvas.draw()
    bb = leg.get_window_extent(ax.figure.canvas.get_renderer())
    n = 0
    for ln in ax.get_lines():
        xy = np.asarray(ln.get_xydata(), float)
        if xy.size == 0:
            continue
        d = ax.transData.transform(xy[np.isfinite(xy).all(axis=1)])
        n += int(np.sum((d[:, 0] > bb.x0) & (d[:, 0] < bb.x1) & (d[:, 1] > bb.y0) & (d[:, 1] < bb.y1)))
    return n


def edge_clipped(png, rim=2, dark=160):
    """Dark pixels on the outer `rim` pixels of a saved PNG (artwork or text cut at the canvas edge)."""
    from PIL import Image
    im = np.asarray(Image.open(png).convert("L"))
    edge = np.concatenate([im[:rim].ravel(), im[-rim:].ravel(), im[:, :rim].ravel(), im[:, -rim:].ravel()])
    return int((edge < dark).sum())


def winkel_x_cm(lon, lat, width_cm):
    """x (cm) of (lon, lat) on a Winkel-Tripel map of total width `width_cm` centred on lon 0 (pygmt R0/<w>c)."""
    def wx(lo, la):
        l, p = np.radians(lo), np.radians(la)
        p1 = np.arccos(2.0 / np.pi)
        a = np.arccos(np.cos(p) * np.cos(l / 2.0))
        return 0.5 * (l * np.cos(p1) + 2.0 * np.cos(p) * np.sin(l / 2.0) / np.sinc(a / np.pi))
    return float(width_cm / 2.0 * (1.0 + wx(lon, lat) / wx(180.0, 0.0)))


def label_left_cm(width_cm, lon, lat, dx_cm, xshift_cm):
    """Page x (cm) of the left end of a left-justified label anchored at (lon, lat) and offset dx_cm, on a panel whose
    left edge is xshift_cm from the page origin."""
    return xshift_cm + winkel_x_cm(lon, lat, width_cm) + dx_cm


def _self_test():
    import os
    import tempfile
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from PIL import Image
    fig, ax = plt.subplots(figsize=(3, 2))
    ax.text(0.1, 0.5, "AAAA"); ax.text(0.12, 0.5, "BBBB")
    assert text_collisions(fig), "overlapping text must be caught"
    plt.close(fig)
    # a tick label outside the view is not drawn and must not count
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(3, 3)); ax1.set_ylim(-500, 400); ax1.set_yticks([-600, -400, 0, 400])
    ax2.set_ylim(0, 10)
    assert not [p for p in text_collisions(fig) if "600" in p[0] or "600" in p[1]], "an out-of-view tick label was counted"
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(4, 3))
    # a dense line along the top edge, fixed limits: it passes under an upper-right legend whatever the font metrics
    ax.plot(np.linspace(0.0, 10.0, 400), np.full(400, 9.5), label="line"); ax.set_xlim(0, 10); ax.set_ylim(0, 10)
    leg = ax.legend(loc="upper right")
    assert legend_over_data(ax, leg) > 0, "a line under the legend must be caught"
    leg.remove(); leg = ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0))
    assert legend_over_data(ax, leg) == 0, "a legend outside the axes covers no data"
    plt.close(fig)
    d = tempfile.mkdtemp()
    a = np.full((50, 80), 255, np.uint8); a[20:30, 30:50] = 0
    Image.fromarray(a).save(os.path.join(d, "ok.png"))
    b = a.copy(); b[10:20, 0:5] = 0
    Image.fromarray(b).save(os.path.join(d, "cut.png"))
    assert edge_clipped(os.path.join(d, "ok.png")) == 0 and edge_clipped(os.path.join(d, "cut.png")) > 0
    assert abs(winkel_x_cm(0.0, 0.0, 10.0) - 5.0) < 1e-9 and abs(winkel_x_cm(-180.0, 0.0, 10.0)) < 1e-9
    assert label_left_cm(7.0, -170.0, 65.0, -1.8, 0.0) < -0.4, "a label 1.8 cm left of 170 W on a 7 cm panel is cut off"
    assert label_left_cm(10.0, -170.0, 65.0, -1.8, 0.0) > 0.0, "the same label on a 10 cm panel is not"
    assert label_left_cm(7.0, -170.0, 65.0, -1.8, 1.2) > 0.5
    print("[figcheck] self-test passed")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--self-test", action="store_true")
    if ap.parse_args().self_test:
        _self_test()
    sys.exit(0)
