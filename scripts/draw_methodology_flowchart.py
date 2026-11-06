"""
=============================================================================
draw_methodology_flowchart.py  —  Figure 1 of the Earth-Science Reviews paper draft
=============================================================================

Renders the methodology flow chart as a clean publication PNG using
matplotlib patches + arrows.  Three colour-coded swim-lanes:

  inputs  (grey)    →    workflow stages (blue)    →    outputs (green)

USAGE
    cd <project>/scripts
    python draw_methodology_flowchart.py
=============================================================================
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from paths_scotese import FIGURES_DIR
import assimilate_scotese as A      # the parameters printed in the figure are read from the workflow itself

# A raster preview of the workflow diagram, written to FIGURES_DIR. The paper's Fig. 2 is the vector diagram drawn by
# esr_revision_flowchart.py.
FIG_DIR = FIGURES_DIR
FIG_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
assert A.DROP_CRATON_SAMPLES and A.DROP_SEA_SITED_SAMPLES and A.DELTA_SMOOTH_SIGMA_FINAL == 0 and A.RESID_REACH_KM[1] > 0, \
    "the flowchart describes the product of record; its settings are not the ones in assimilate_scotese.py"

# ---------------------------------------------------------------------------
# Layout — coordinates in 0..1 figure units
# ---------------------------------------------------------------------------
COLORS = {
    "input":  "#dcdcdc",   # grey
    "stage":  "#cfd9e6",   # light blue
    "iter":   "#ffe4b5",   # light orange (the per-slice loop body)
    "output": "#cfe6cf",   # light green
    "edge":   "#444444",
}

def add_box(ax, x, y, w, h, text, color, fontsize=9, weight="normal"):
    box = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.012,rounding_size=0.012",
        facecolor=color, edgecolor="black", linewidth=0.9,
    )
    ax.add_patch(box)
    ax.text(x + w/2, y + h/2, text, ha="center", va="center",
            fontsize=fontsize, weight=weight, wrap=True)


def arrow(ax, x0, y0, x1, y1, style="->", lw=1.3, color="#444"):
    a = FancyArrowPatch(
        (x0, y0), (x1, y1),
        arrowstyle=style, mutation_scale=14, lw=lw, color=color,
    )
    ax.add_patch(a)


fig, ax = plt.subplots(figsize=(13, 16))
ax.set_xlim(0, 1); ax.set_ylim(0, 1)
ax.set_axis_off()

# ---- Lane headers ------------------------------------------------------
ax.text(0.18, 0.945, "INPUTS",          ha="center", fontsize=11, weight="bold", color="#555")
ax.text(0.50, 0.945, "ASSIMILATION",    ha="center", fontsize=11, weight="bold", color="#555")
ax.text(0.82, 0.945, "OUTPUTS",         ha="center", fontsize=11, weight="bold", color="#555")

# ---- INPUT boxes (left column) ----------------------------------------
add_box(ax, 0.02, 0.86, 0.32, 0.05,
        "Kinematic paleo-topography\n"
        "(Scotese & Wright 2018 PaleoDEMs\n"
        "or Merdith et al. 2024)",
        COLORS["input"])

add_box(ax, 0.02, 0.76, 0.32, 0.06,
        "Geochemical sample database\n"
        "(~100k igneous rocks; major + trace elements,\n"
        "present-day Lat/Lon, sample age, Tecto_Prov)",
        COLORS["input"])

add_box(ax, 0.02, 0.66, 0.32, 0.06,
        "Plate model (rotation file, present-day static\n"
        "continental polygons; plate-boundary topologies\n"
        "for 0–100 Ma only, used for the arc province)",
        COLORS["input"])

# ---- WORKFLOW boxes (middle column) -----------------------------------
# One-time preprocessing
add_box(ax, 0.36, 0.86, 0.30, 0.05,
        "Crustal-thickness mohometry\n(Zhou et al. 2025, ML on 66 features)",
        COLORS["stage"], fontsize=9)

add_box(ax, 0.36, 0.78, 0.30, 0.05,
        "Airy isostasy → 6-model ensemble\nper-sample elevation z ± σ",
        COLORS["stage"], fontsize=9)

add_box(ax, 0.36, 0.70, 0.30, 0.05,
        "Plate-ID assignment\n"
        "(point-in-polygon against\nstatic continental polygons, anchor 000)",
        COLORS["stage"], fontsize=9)

add_box(ax, 0.36, 0.615, 0.30, 0.05,
        "Drop samples not in any polygon (~2.5%)\n"
        "and samples inside craton outlines\n(Shirmard et al. 2025)",
        COLORS["stage"], fontsize=9)

# Per-slice loop banner
add_box(ax, 0.30, 0.555, 0.42, 0.04,
        "For each target age slice t ∈ {0, 5, 10, …}:",
        COLORS["iter"], fontsize=10, weight="bold")

# Per-slice steps
add_box(ax, 0.36, 0.49, 0.30, 0.04,
        "1.  Sample window |Age − t| ≤ Δt(t),\ntriangular weight in |Age − t|",
        COLORS["stage"], fontsize=9)

add_box(ax, 0.36, 0.43, 0.30, 0.04,
        "2.  Reconstruct sample positions to age t\n(via plate ID + rotation model)",
        COLORS["stage"], fontsize=9)

add_box(ax, 0.36, 0.36, 0.30, 0.04,
        "3.  Drop samples landing outside\nthe per-slice continent mask",
        COLORS["stage"], fontsize=9)

add_box(ax, 0.36, 0.30, 0.30, 0.04,
        "4.  Decluster (1° × Tecto_Prov bins, weighted\nmedian + Kish); drop bins on the prior's sea",
        COLORS["stage"], fontsize=9)

add_box(ax, 0.36, 0.24, 0.30, 0.04,
        "5.  Build tectonic-province grid\n(1-cell speckle removed)",
        COLORS["stage"], fontsize=9)

add_box(ax, 0.36, 0.18, 0.30, 0.04,
        f"6.  Province-wise CDF rescaling (shrunk N_eff {A.N_MIN_P}→{A.N_FULL_P}),\n"
        f"zero beyond {A.LOCAL_SUPPORT_KM:.0f} km of the nearest sample,\n"
        f"province-edge taper, cross-boundary σ {A.CROSS_BOUNDARY_SIGMA:g} cells",
        COLORS["stage"], fontsize=9)

add_box(ax, 0.36, 0.12, 0.30, 0.04,
        f"7.  Smoothed-residual kernel ({A.RESID_LS_KM:.0f} km Gaussian),\n"
        f"whole to {A.RESID_REACH_KM[0]:.0f} km of the nearest sample, zero at {A.RESID_REACH_KM[1]:.0f} km",
        COLORS["stage"], fontsize=9)

add_box(ax, 0.36, 0.06, 0.30, 0.04,
        f"8.  Continent-masked Δz smoothing (σ {A.DELTA_SMOOTH_SIGMA:g} cells)\n"
        f"+ cap |Δz| ≤ {A.DZ_MAX/1000:g} km; correction ramped to 0 at 0 Ma\n"
        f"over the youngest {A.YOUNG_TAPER_T_FULL:g} Myr",
        COLORS["stage"], fontsize=9)

# ---- OUTPUT boxes (right column) --------------------------------------
add_box(ax, 0.68, 0.86, 0.30, 0.05,
        "Corrected per-Ma NetCDF\n(M_orig, M_corrected, Δz, provinces,\n"
        "continent_mask, n_eff per cell)",
        COLORS["output"], fontsize=9)

add_box(ax, 0.68, 0.76, 0.30, 0.05,
        "Per-slice statistics CSV\n"
        "(bias before/after, RMS, hypsometric pcts,\n"
        " offshore-drop rate, …)",
        COLORS["output"], fontsize=9)

add_box(ax, 0.68, 0.66, 0.30, 0.04,
        "Province-stratified Δz CSV",
        COLORS["output"], fontsize=9)

add_box(ax, 0.68, 0.60, 0.30, 0.04,
        "Era-binned summary markdown",
        COLORS["output"], fontsize=9)

add_box(ax, 0.68, 0.51, 0.30, 0.06,
        "Diagnostic figures\n"
        "(time-series of bias, RMS, p99 land,\n"
        "hypsometric heatmaps, per-province Δz)",
        COLORS["output"], fontsize=9)

add_box(ax, 0.68, 0.39, 0.30, 0.08,
        "Four video products per dataset\n"
        "(Winkel-Tripel, fixed cpt ±5 km, end arrows):\n"
        "1. original elevation\n"
        "2. corrected elevation\n"
        "3. Δz (corrected – original)\n"
        "4. corrected + geochem control points",
        COLORS["output"], fontsize=9)

add_box(ax, 0.68, 0.26, 0.30, 0.06,
        "Application products:\n"
        "supermountains index, area-weighted\n"
        "hypsometric curves, etc.",
        COLORS["output"], fontsize=9)

# ---- Arrows — inputs → preprocessing -----------------------------------
# geochem → mohometry
arrow(ax, 0.34, 0.79, 0.36, 0.885)
arrow(ax, 0.34, 0.79, 0.36, 0.805)
# plate model → plate-ID assignment + reconstruction
arrow(ax, 0.34, 0.69, 0.36, 0.725)

# preprocessing chain
arrow(ax, 0.51, 0.86, 0.51, 0.83)   # mohometry → isostasy
arrow(ax, 0.51, 0.78, 0.51, 0.75)   # isostasy → plate-ID
arrow(ax, 0.51, 0.70, 0.51, 0.665)  # plate-ID → drop

# preprocessing → per-slice loop
arrow(ax, 0.51, 0.615, 0.51, 0.595)

# per-slice loop banner → step 1
arrow(ax, 0.51, 0.555, 0.51, 0.53)

# step chain (within per-slice)
for y0, y1 in [(0.49, 0.47),
               (0.43, 0.40),
               (0.36, 0.34),
               (0.30, 0.28),
               (0.24, 0.22),
               (0.18, 0.16),
               (0.12, 0.10)]:
    arrow(ax, 0.51, y0, 0.51, y1)

# Per-slice kinematic input → step 1
arrow(ax, 0.34, 0.88, 0.36, 0.51)

# Final step → outputs
arrow(ax, 0.66, 0.08, 0.78, 0.41, color="#1a6f1a")    # to videos
arrow(ax, 0.66, 0.08, 0.78, 0.51, color="#1a6f1a")    # to diagnostics
arrow(ax, 0.66, 0.08, 0.68, 0.62, color="#1a6f1a")    # to era summary
arrow(ax, 0.66, 0.08, 0.68, 0.68, color="#1a6f1a")    # to province CSV
arrow(ax, 0.66, 0.08, 0.68, 0.78, color="#1a6f1a")    # to per-slice CSV
arrow(ax, 0.66, 0.08, 0.68, 0.88, color="#1a6f1a")    # to corrected NetCDF
arrow(ax, 0.66, 0.08, 0.78, 0.27, color="#1a6f1a")    # to application

# Loop-back arrow indicating "for each slice"
arrow(ax, 0.30, 0.575, 0.30, 0.08,
      style="<-", color="#aa6e10", lw=1.0)
ax.text(0.27, 0.30, "loop:\n0 → 540 Ma\n(every 5 Myr)",
        rotation=90, fontsize=8.5, color="#aa6e10",
        ha="right", va="center")

# ---- footer ------------------------------------------------------------
ax.text(0.5, 0.02,
        "Parameters of the product:  "
        f"Δt(t) = 5/10/20 Myr   |   N_min,p = {A.N_MIN_P}, N_full,p = {A.N_FULL_P}   |   "
        f"support {A.LOCAL_SUPPORT_KM:.0f} km   |   ℓ_resid = {A.RESID_LS_KM:.0f} km, reach "
        f"{A.RESID_REACH_KM[0]:.0f}–{A.RESID_REACH_KM[1]:.0f} km   |   "
        f"|Δz| ≤ {A.DZ_MAX/1000:g} km   |   σ_smooth = {A.DELTA_SMOOTH_SIGMA:g} cells",
        ha="center", fontsize=8.5, style="italic", color="#444")

plt.tight_layout()
out_png  = FIG_DIR    / "methodology_flowchart_v1.png"
plt.savefig(out_png, dpi=180, bbox_inches="tight")
print(f"wrote {out_png}")
