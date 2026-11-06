"""
=============================================================================
sensitivity_refinements.py  —  Ablation study of the four methodological
refinements added on top of the baseline assimilation
=============================================================================

For each of N representative ages, the script runs `assimilate_one` under
five configurations:

  1. Default                 — all four refinements ON
  2. No temporal kernel       — uniform weighting inside Δt(t)
  3. Hard N=30 threshold      — no smooth shrinkage
  4. No province pooling      — each province uses only its own samples
  5. All refinements OFF      — legacy behaviour (uniform + hard + no pool)

For each (age × configuration) it records the headline diagnostics
(bias before/after, RMS before/after, p99 land elevation, Δz RMS,
declustered sample count) and produces a comparison figure suitable
for the paper supplementary.

OUTPUTS
    OUTPUT_DIR/tables/sensitivity_refinements.csv
    FIGURES_DIR/sensitivity_refinements.{png,pdf} (Fig. S4; build_all_figures.sh gives it its number); panels
    lettered, no title,
    refused if any text overlaps (figure_checks.text_collisions)

USAGE
    cd <project>/scripts
    python sensitivity_refinements.py                  # six default ages
    python sensitivity_refinements.py --ages 50 100 200 400 540
    python sensitivity_refinements.py --plot-only      # redraw from the CSV

Runtime: ~5 ages × 5 configs × ~2 s/slice ≈ 1 minute.
=============================================================================
"""
from __future__ import annotations
import argparse, sys, copy
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import assimilate_scotese as A
from paths_scotese import FIGURES_DIR, OUTPUT_DIR
from figure_checks import text_collisions
PROJ_ROOT = HERE.parent
FIG_DIR = FIGURES_DIR
FIG_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Configuration matrix
# ---------------------------------------------------------------------------
CONFIGS = [
    ("Product (all refinements on)",
     dict(TEMPORAL_KERNEL="triangular", N_MIN_P=3,  N_FULL_P=30, PROVINCE_POOLING=True),
     "#1f77b4"),
    ("No temporal kernel (uniform weights)",
     dict(TEMPORAL_KERNEL="uniform",    N_MIN_P=3,  N_FULL_P=30, PROVINCE_POOLING=True),
     "#ff7f0e"),
    ("Hard threshold N = 30 (no shrinkage)",
     dict(TEMPORAL_KERNEL="triangular", N_MIN_P=30, N_FULL_P=30, PROVINCE_POOLING=True),
     "#2ca02c"),
    ("No province pooling",
     dict(TEMPORAL_KERNEL="triangular", N_MIN_P=3,  N_FULL_P=30, PROVINCE_POOLING=False),
     "#d62728"),
    ("All four off",
     dict(TEMPORAL_KERNEL="uniform",    N_MIN_P=30, N_FULL_P=30, PROVINCE_POOLING=False),
     "#7f7f7f"),
]

DEFAULT_AGES = [50, 100, 200, 300, 400, 500]


def _set_config(cfg: dict):
    """Apply a configuration to assimilate_scotese module globals."""
    for key, value in cfg.items():
        setattr(A, key, value)


def _snapshot_config():
    """Capture current values of the four toggleable parameters."""
    return {k: getattr(A, k) for k in
            ("TEMPORAL_KERNEL", "N_MIN_P", "N_FULL_P", "PROVINCE_POOLING")}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--ages", type=int, nargs="+", default=DEFAULT_AGES,
                   help="ages (Ma) to test (default: 50 100 200 300 400 500)")
    p.add_argument("--plot-only", action="store_true", help="redraw the figure from the CSV")
    args = p.parse_args()
    out_csv = OUTPUT_DIR / "tables" / "sensitivity_refinements.csv"
    if args.plot_only:
        plot(pd.read_csv(out_csv))
        return

    # Stash original config so we can restore it
    original_cfg = _snapshot_config()
    print(f"Original config (will be restored at end):\n  {original_cfg}\n")

    rows = []
    for label, cfg, color in CONFIGS:
        _set_config(cfg)
        print(f"--- config: {label} ---")
        for t in args.ages:
            summary, *_ = A.assimilate_one(float(t), save_nc=False)
            summary["config"] = label
            summary["color"] = color
            rows.append(summary)
            print(f"   t={t:3d} Ma  bias {summary['bias_before_m']:6.0f} → "
                  f"{summary['bias_after_m']:6.0f},  "
                  f"RMS {summary['rms_before_m']:6.0f} → {summary['rms_after_m']:6.0f},  "
                  f"p99 {summary['p99_after_m']:6.0f},  Δrms {summary['delta_rms_m']:6.0f},  "
                  f"n_dec {summary['n_decluster']}")
        print()

    # Restore original config so subsequent imports of assimilate_scotese
    # don't see the last test setting
    _set_config(original_cfg)

    # Every configuration must change something: a configuration identical to Default
    # means the setting did not reach assimilate_scotese (for example a module constant
    # bound as a default argument at import), and is an error, not a result.
    cols = ["bias_after_m", "rms_after_m", "p99_after_m", "delta_rms_m"]
    ref = pd.DataFrame([r for r in rows if r["config"] == CONFIGS[0][0]]).set_index("t_Ma")[cols]
    for label, _, _ in CONFIGS[1:]:
        cur = pd.DataFrame([r for r in rows if r["config"] == label]).set_index("t_Ma")[cols]
        dmax = float(np.nanmax(np.abs(cur.values - ref.loc[cur.index].values)))
        print(f"  {label}: largest difference from Default over the four metrics {dmax:.1f} m")
        if dmax == 0.0:
            sys.exit(f"'{label}' is identical to Default at every age: the toggle did not reach "
                     "assimilate_scotese. Nothing written.")

    df = pd.DataFrame(rows)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    print(f"wrote {out_csv}")
    plot(df)


def plot(df):
    metrics = [
        ("bias_after_m", "sample bias after correction (m)"),
        ("rms_after_m", "sample RMS after correction (m)"),
        ("p99_after_m", "99th percentile of land elevation (m)"),
        ("delta_rms_m", "RMS of \u0394z over the continents (m)"),
    ]
    plt.rcParams.update({"font.size": 10})
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    for k, (ax, (col, lab)) in enumerate(zip(axes.flat, metrics)):
        for label, _, color in CONFIGS:
            sub = df[df["config"] == label].sort_values("t_Ma")
            first = label == CONFIGS[0][0]
            ls = "-" if first else (":" if label == CONFIGS[-1][0] else "--")
            ax.plot(sub["t_Ma"], sub[col], "o" + ls, color=color, lw=2.2 if first else 1.4, ms=5, label=label)
        ax.set_xlabel("age (Ma)")
        ax.set_ylabel(lab)
        ax.invert_xaxis()
        ax.grid(True, alpha=0.3, lw=0.5)
        ax.text(-0.13, 1.03, "abcd"[k], transform=ax.transAxes, fontsize=13, fontweight="bold", va="bottom")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.legend(handles, labels, loc="upper center", ncol=3, fontsize=9.5, bbox_to_anchor=(0.5, 1.0), frameon=False)
    bad = text_collisions(fig)
    if bad:
        raise SystemExit(f"[sensitivity] overlapping text {bad[:6]}")
    for ext in ("png", "pdf"):
        fig.savefig(FIG_DIR / f"sensitivity_refinements.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {FIG_DIR / 'sensitivity_refinements'}.png/.pdf (text-collision check passed)")


if __name__ == "__main__":
    main()
