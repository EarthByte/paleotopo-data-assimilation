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
    data/corrected/sensitivity_refinements.csv
    Figures/Fig10_sensitivity_refinements.png

USAGE
    cd <project>/scripts_Scotese
    python sensitivity_refinements.py                  # six default ages
    python sensitivity_refinements.py --ages 50 100 200 400 540

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
from paths_scotese import CORRECTED_DIR
PROJ_ROOT = HERE.parent
FIG_DIR = PROJ_ROOT / "Figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Configuration matrix
# ---------------------------------------------------------------------------
CONFIGS = [
    ("Default (all refinements on)",
     dict(TEMPORAL_KERNEL="triangular", N_MIN_P=3,  N_FULL_P=30, PROVINCE_POOLING=True),
     "#1f77b4"),
    ("(1) No temporal kernel",
     dict(TEMPORAL_KERNEL="uniform",    N_MIN_P=3,  N_FULL_P=30, PROVINCE_POOLING=True),
     "#ff7f0e"),
    ("(2,3) Hard N=30 threshold (no shrinkage)",
     dict(TEMPORAL_KERNEL="triangular", N_MIN_P=30, N_FULL_P=30, PROVINCE_POOLING=True),
     "#2ca02c"),
    ("(4) No province pooling",
     dict(TEMPORAL_KERNEL="triangular", N_MIN_P=3,  N_FULL_P=30, PROVINCE_POOLING=False),
     "#d62728"),
    ("All refinements off (legacy)",
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
    args = p.parse_args()

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

    # ----------------------------------------------------------------------
    # Save CSV
    # ----------------------------------------------------------------------
    df = pd.DataFrame(rows)
    out_csv = CORRECTED_DIR / "sensitivity_refinements.csv"
    df.to_csv(out_csv, index=False)
    print(f"wrote {out_csv}")

    # ----------------------------------------------------------------------
    # Plot
    # ----------------------------------------------------------------------
    metrics = [
        ("bias_after_m",  "Sample bias after correction (m)"),
        ("rms_after_m",   "Sample-residual RMS after correction (m)"),
        ("p99_after_m",   "p99 land elevation after correction (m)"),
        ("delta_rms_m",   "Δz RMS over continent (m)"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(13, 8.5))
    for ax, (metric_col, metric_label) in zip(axes.flat, metrics):
        for label, cfg, color in CONFIGS:
            sub = df[df["config"] == label].sort_values("t_Ma")
            ls = "-" if "Default" in label else (":" if "legacy" in label else "--")
            lw = 2.2 if "Default" in label else 1.4
            ax.plot(sub["t_Ma"], sub[metric_col], "o" + ls,
                    color=color, lw=lw, ms=5, label=label)
        ax.set_xlabel("Age (Ma)")
        ax.set_ylabel(metric_label, fontsize=10)
        ax.invert_xaxis()
        ax.grid(True, alpha=0.3, lw=0.5)
    # Single legend at the top
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, fontsize=9,
               bbox_to_anchor=(0.5, 1.02), frameon=True)
    plt.suptitle("Ablation of the four methodological refinements",
                 fontsize=13, y=1.06)
    plt.tight_layout()

    out_png = FIG_DIR / "Fig10_sensitivity_refinements.png"
    plt.savefig(out_png, dpi=180, bbox_inches="tight")
    print(f"wrote {out_png}")


if __name__ == "__main__":
    main()
