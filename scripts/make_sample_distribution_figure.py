"""
=============================================================================
make_sample_distribution_figure.py  —  Figure 2 of the Earth-Science Reviews paper
=============================================================================

Four-panel figure characterising the geochemical sample database used in
the assimilation:

  (a) Global map of present-day sample positions, coloured by Age_Ma.
  (b) Histogram of sample count per 10 Myr age bin (Phanerozoic, 0–540 Ma).
  (c) Schematic of the per-sample elevation calculation:
      whole-rock geochemistry → ML mohometry → 6-model Airy isostasy
      ensemble → ensemble mean + σ.
  (d) Per-province declustered sample counts, stratified by era.

OUTPUT
    paper/Scotese/Fig02_sample_distribution.png

USAGE
    cd <project>/scripts_Scotese
    python make_sample_distribution_figure.py
=============================================================================
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from paths_scotese import CSV_PATH

PROJECT_ROOT = HERE.parent
FIG_DIR = PROJECT_ROOT / "paper" / "Scotese"
FIG_DIR.mkdir(parents=True, exist_ok=True)


def era_of(t):
    if t < 66: return "Cenozoic"
    if t < 252: return "Mesozoic"
    if t < 540: return "Paleozoic"
    return "older"


def main():
    print("Loading geochem CSV …")
    df = pd.read_csv(CSV_PATH, low_memory=False)
    df = df.dropna(subset=["Lat", "Lon", "Age_Ma"]).copy()
    # Restrict to the same Phanerozoic window the pipeline uses
    df = df[(df["Age_Ma"] >= 0) & (df["Age_Ma"] <= 540)].reset_index(drop=True)
    print(f"  {len(df)} samples in 0–540 Ma")

    fig = plt.figure(figsize=(13, 9.5))

    # ----------------------- (a) global map -----------------------
    ax_a = fig.add_subplot(2, 2, 1, projection=ccrs.Robinson())
    ax_a.set_global()
    ax_a.coastlines(lw=0.4, color="grey")
    ax_a.gridlines(lw=0.3, color="grey", alpha=0.3)
    sc = ax_a.scatter(df["Lon"], df["Lat"], c=df["Age_Ma"],
                      s=5, alpha=0.45, cmap="plasma_r",
                      vmin=0, vmax=540,
                      transform=ccrs.PlateCarree(), edgecolor="none")
    cbar = fig.colorbar(sc, ax=ax_a, orientation="horizontal",
                        fraction=0.05, pad=0.06)
    cbar.set_label("sample age (Ma)", fontsize=9)
    ax_a.set_title(f"(a) Present-day positions of {len(df):,} samples (0–540 Ma)",
                   fontsize=11)

    # ----------------------- (b) age histogram -----------------------
    ax_b = fig.add_subplot(2, 2, 2)
    bins = np.arange(0, 545, 10)
    n, edges, _ = ax_b.hist(df["Age_Ma"], bins=bins,
                            color="#3a6c9e", edgecolor="white", lw=0.5)
    ax_b.set_xlabel("sample age (Ma)"); ax_b.set_ylabel("samples per 10 Myr bin")
    ax_b.set_title(f"(b) Sample density vs age", fontsize=11)
    ax_b.invert_xaxis()
    ax_b.set_yscale("log"); ax_b.grid(True, alpha=0.3)
    # annotate era boundaries
    for boundary, name in [(66, "K-Pg"), (252, "P-Tr")]:
        ax_b.axvline(boundary, color="#444", lw=0.5, ls="--")
        ax_b.text(boundary, 1.02, name, transform=ax_b.get_xaxis_transform(),
                  ha="center", va="bottom", fontsize=8, color="#666")

    # ----------------------- (c) elevation calc schematic -----------------------
    ax_c = fig.add_subplot(2, 2, 3); ax_c.set_axis_off()
    def box(ax, x, y, w, h, txt, color="#cfd9e6", fs=9, bold=False):
        from matplotlib.patches import FancyBboxPatch
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012",
                                    facecolor=color, edgecolor="black", lw=0.8))
        ax.text(x + w/2, y + h/2, txt, ha="center", va="center", fontsize=fs,
                weight="bold" if bold else "normal")
    def arr(ax, x0, y0, x1, y1):
        from matplotlib.patches import FancyArrowPatch
        ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="->",
                                     mutation_scale=12, lw=1.1, color="#444"))
    ax_c.set_xlim(0, 1); ax_c.set_ylim(0, 1)
    box(ax_c, 0.05, 0.80, 0.90, 0.10,
        "Whole-rock geochemistry\n(66 major + trace features per sample)",
        color="#dcdcdc")
    arr(ax_c, 0.50, 0.80, 0.50, 0.72)
    box(ax_c, 0.05, 0.59, 0.90, 0.13,
        "Machine-learning mohometry\n(random forest; Zhou et al. 2025)\n→ crustal thickness  +  uncertainty envelope",
        color="#cfd9e6")
    arr(ax_c, 0.50, 0.59, 0.50, 0.50)
    box(ax_c, 0.05, 0.32, 0.90, 0.18,
        "Airy isostasy with 6 density models\n(Herzberg 0.08 / 0.23 / 0.38, Brown 2022,\nDavis 2009, Condie 2016)",
        color="#cfd9e6")
    arr(ax_c, 0.50, 0.32, 0.50, 0.24)
    box(ax_c, 0.05, 0.05, 0.90, 0.18,
        "Per-sample observed elevation z ± σ\n"
        "σ² = σ²(ensemble spread) + σ²(thickness envelope) + σ²(quality)\n"
        "  →  inverse-variance weight in the assimilation",
        color="#cfe6cf")
    ax_c.set_title("(c) Per-sample elevation calculation", fontsize=11, pad=14)

    # ----------------------- (d) declustered counts by province × era -----------------------
    ax_d = fig.add_subplot(2, 2, 4)
    df["era"] = df["Age_Ma"].apply(era_of)
    df["prov"] = df.get("Tecto_Prov", pd.Series("Other", index=df.index)).fillna("Other")
    pivot = df.pivot_table(index="prov", columns="era", values="Lat",
                           aggfunc="count", fill_value=0)
    desired_provs = ["Orogen", "Continental Arc", "Continental Margin",
                     "Extended Crust", "Island Arc", "Basin", "Platform",
                     "Shield", "Other"]
    pivot = pivot.reindex([p for p in desired_provs if p in pivot.index])
    eras = ["Cenozoic", "Mesozoic", "Paleozoic"]
    pivot = pivot[[e for e in eras if e in pivot.columns]]
    pivot.plot.barh(ax=ax_d, stacked=True, width=0.72,
                    color=["#fae6c8", "#cfe6cf", "#cfd9e6"][:pivot.shape[1]],
                    edgecolor="black", lw=0.3)
    ax_d.set_xlabel("raw sample count")
    ax_d.set_ylabel("")
    ax_d.invert_yaxis()
    ax_d.grid(True, axis="x", alpha=0.3)
    ax_d.legend(title="Era", fontsize=9, framealpha=0.9)
    ax_d.set_title("(d) Sample count by tectonic province × era", fontsize=11)

    plt.tight_layout()
    out_png = FIG_DIR / "Fig02_sample_distribution.png"
    plt.savefig(out_png, dpi=200, bbox_inches="tight")
    print(f"  wrote {out_png}")


if __name__ == "__main__":
    main()
