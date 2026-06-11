#!/usr/bin/env python3
"""
=============================================================================
make_flooded_fraction_figure.py  —  Continental flooding summary stats:
                                    corrected vs corrected + per-step dyntopo increment
=============================================================================

WHAT THIS DOES
    For each requested age, computes the percentage of CONTINENTAL grid
    cells (continent_mask == True) where elevation is below sea level
    (z < 0 m) in two fields:

        (1) the geochem-corrected Scotese & Wright paleotopography
            `M_corrected` - the main-text result
        (2) the same map with the Young 2022 dyntopo PER-STEP INCREMENT
            relative to present added on top
            `M_combined = M_corrected + [z_dyntopo(t) - z_dyntopo(t - Δt)]`

    Both percentages are plotted vs age on a simple xy line plot, so
    the reader can see at a glance whether composing the per-step dyntopo increment into the corrected map drives
    continental flooding outside
    the bounds of what's geologically reasonable (cf. Kocsis & Scotese,
    ESR, 2021).

    Default: ALL ages with a combined NetCDF in --combined-dir (the
    full 0-300 Ma / 5-Myr Young 2022 GLD428 availability window after
    `add_dyntopo_to_corrected_scotese.py` has run, ~61 points).  The
    plot is dense and shows the temporal evolution of the dyntopo
    contribution to continental flooding.  Override with --ages to
    evaluate a sparse subset if you only want anchor points.

OUTPUT
    paths_scotese.OUTPUT_DIR / "SW_flooded_fraction_corrected_vs_dyntopo.png"
    paths_scotese.OUTPUT_DIR / "SW_flooded_fraction_corrected_vs_dyntopo.pdf"
    paths_scotese.OUTPUT_DIR / "SW_flooded_fraction_corrected_vs_dyntopo.csv"
        (the same numbers as a small csv, for traceability)

USAGE
    cd <project>/scripts_Scotese
    # Make sure the combined NetCDFs are present first:
    python add_dyntopo_to_corrected_scotese.py --dyntopo-dir <...>
    python make_flooded_fraction_figure.py
    python make_flooded_fraction_figure.py --ages 410 300 250 150 75 25

OPTIONS
    --ages INT ...              ages to evaluate (default 500 400 300 200 100 50)
    --combined-dir PATH         directory of <age>Ma_corrected_plus_dyntopo_SW.nc
                                (default: PROJECT_ROOT/data/corrected_Scotese_plus_dyntopo)
    --no-pdf                    skip PDF export
    --no-png                    skip PNG export
    --no-csv                    skip CSV export
    --dpi INT                   PNG dpi (default 200)
    --figsize FLOAT FLOAT       (width, height) in inches (default 7 4.5)

DEPENDENCIES
    matplotlib, netCDF4, numpy
=============================================================================
"""
from __future__ import annotations
import argparse, sys, csv
from pathlib import Path
import numpy as np
import netCDF4 as nc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from paths_scotese import OUTPUT_DIR, PROJECT_ROOT


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
import re
COMBINED_FNAME_FMT = "{age}Ma_corrected_plus_dyntopo_diff_young_SW.nc"
COMBINED_FNAME_RE  = re.compile(r"^(\d+)Ma_corrected_plus_dyntopo_diff_young_SW\.nc$")
DEFAULT_COMBINED_DIR = PROJECT_ROOT / "data" / "corrected_Scotese_plus_dyntopo_diff_young"


def auto_discover_ages(combined_dir: Path) -> list[int]:
    """Return sorted list of ages with a combined NetCDF in `combined_dir`."""
    if not combined_dir.is_dir():
        return []
    return sorted(
        int(COMBINED_FNAME_RE.match(p.name).group(1))
        for p in combined_dir.iterdir()
        if COMBINED_FNAME_RE.match(p.name)
    )


# ---------------------------------------------------------------------------
# Flooded-fraction calculation
# ---------------------------------------------------------------------------
def flooded_fraction(M: np.ndarray, cont_mask: np.ndarray) -> float:
    """Return percent of continental cells with elevation < 0.

    A "flooded continental cell" is any (lat, lon) where
        continent_mask is True AND M < 0.
    Result is in percent of total continental cells (NaN-safe: NaN cells
    are excluded from both numerator and denominator).
    """
    cont_mask = cont_mask.astype(bool)
    cells_continental = cont_mask & np.isfinite(M)
    n_total = int(cells_continental.sum())
    if n_total == 0:
        return float("nan")
    n_flooded = int((cells_continental & (M < 0)).sum())
    return 100.0 * n_flooded / n_total


def load_age(combined_dir: Path, age: int):
    """Read M_corrected, M_combined, and the continent_mask for `age`.

    M_corrected and M_combined come from the build_dyntopo_diff_correction.py
    output NetCDF.  The continent mask is NOT carried in those NetCDFs;
    we read it from the corresponding corrected-S&W NetCDF instead.
    """
    f = combined_dir / COMBINED_FNAME_FMT.format(age=int(age))
    if not f.exists():
        raise FileNotFoundError(
            f"{f} not found - run `python build_dyntopo_diff_correction.py "
            f"--source young --ages {int(age)}` first")
    with nc.Dataset(f) as d:
        Mc    = np.asarray(d["M_corrected"][:], dtype=float)
        Mcomb = np.asarray(d["M_combined"][:], dtype=float)
    corrected_sw = (PROJECT_ROOT / "data" / "corrected_Scotese"
                    / f"{int(age)}Ma_corrected_SW.nc")
    if corrected_sw.exists():
        with nc.Dataset(corrected_sw) as d:
            cont = np.asarray(d["continent_mask"][:], dtype=bool)
    else:
        cont = (Mc >= 0)
    return Mc, Mcomb, cont


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
def make_figure(ages: list[int], combined_dir: Path,
                figsize=(7.0, 4.5),
                write_png=True, write_pdf=True, write_csv=True,
                dpi=200):
    print(f"computing flooded-fraction summary for {len(ages)} ages: {ages}")
    rows = []
    for a in sorted(ages, reverse=True):     # oldest -> youngest for the plot
        Mc, Mcomb, cont = load_age(combined_dir, a)
        f_corr = flooded_fraction(Mc,    cont)
        f_comb = flooded_fraction(Mcomb, cont)
        print(f"  {a:4d} Ma:  corrected = {f_corr:5.1f}%   "
              f"corrected+dyntopo = {f_comb:5.1f}%   "
              f"Delta = {f_comb - f_corr:+5.1f} pp")
        rows.append((a, f_corr, f_comb))

    ages_arr  = np.array([r[0] for r in rows])
    f_corr_arr = np.array([r[1] for r in rows])
    f_comb_arr = np.array([r[2] for r in rows])

    # Dense (>= 12 points) vs sparse defaults: shrink markers and drop
    # explicit per-age xticks when we have a full sweep, so the line
    # plot reads as a continuous time series rather than a scatter.
    dense = len(ages) >= 12
    ms = 4 if dense else 8
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(ages_arr, f_corr_arr,
            marker="o", color="#1f77b4", lw=1.8, ms=ms,
            label="Corrected (this study)")
    ax.plot(ages_arr, f_comb_arr,
            marker="s", color="#d62728", lw=1.8, ms=ms, ls="--",
            label="Corrected + Young 2022 dyntopo per-step increment")
    # Plot direction: oldest on the left, youngest on the right
    ax.invert_xaxis()
    ax.set_xlabel("Age (Ma)", fontsize=12)
    ax.set_ylabel("Continental cells below sea level (%)", fontsize=12)
    ax.grid(True, alpha=0.3, linestyle=":")
    ax.legend(loc="best", frameon=False, fontsize=10)
    ax.tick_params(labelsize=10)
    if dense:
        # auto-tick every 50 Ma when running the full sweep
        ax.set_xticks(np.arange(0, int(max(ages)) + 50, 50))
        ax.set_xlim(max(ages) + 10, min(ages) - 10)
    else:
        # sparse anchor mode: explicit ticks at each anchor age
        ax.set_xticks(sorted(ages, reverse=True))
        ax.set_xlim(max(ages) + 25, min(ages) - 25)
    # Set sensible y-limit floor at 0 and pad the top a bit
    ymax = max(np.nanmax(f_corr_arr), np.nanmax(f_comb_arr))
    ax.set_ylim(0, max(10, ymax * 1.15))
    fig.tight_layout()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    out_base = OUTPUT_DIR / "SW_flooded_fraction_corrected_vs_dyntopo"
    if write_png:
        out_png = out_base.with_suffix(".png")
        fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
        print(f"  wrote {out_png}")
    if write_pdf:
        out_pdf = out_base.with_suffix(".pdf")
        fig.savefig(out_pdf, bbox_inches="tight")
        print(f"  wrote {out_pdf}")
    plt.close(fig)
    if write_csv:
        out_csv = out_base.with_suffix(".csv")
        with open(out_csv, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["age_ma",
                        "pct_continent_flooded_corrected",
                        "pct_continent_flooded_corrected_plus_dyntopo_per_step_increment",
                        "delta_percentage_points"])
            for a, fc, fcb in rows:
                w.writerow([a, f"{fc:.2f}", f"{fcb:.2f}", f"{fcb - fc:+.2f}"])
        print(f"  wrote {out_csv}")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ages", type=int, nargs="+", default=None,
                   help="ages to evaluate. Default: ALL ages with a "
                        "combined NetCDF in --combined-dir (full 0-300 Ma / "
                        "5-Myr sweep). Pass an explicit list to evaluate "
                        "anchor points only, e.g. `--ages 300 200 100 50 0`.")
    p.add_argument("--combined-dir", type=Path, default=DEFAULT_COMBINED_DIR,
                   help="directory of <age>Ma_corrected_plus_dyntopo_SW.nc files")
    p.add_argument("--no-png", action="store_true")
    p.add_argument("--no-pdf", action="store_true")
    p.add_argument("--no-csv", action="store_true")
    p.add_argument("--dpi", type=int, default=200)
    p.add_argument("--figsize", type=float, nargs=2, default=(7.0, 4.5),
                   metavar=("WIDTH", "HEIGHT"),
                   help="figure size in inches (default 7 4.5)")
    args = p.parse_args()

    # Auto-discover full sweep when --ages not given
    if args.ages is None:
        args.ages = auto_discover_ages(args.combined_dir)
        if not args.ages:
            raise SystemExit(
                f"ERROR: no <age>Ma_corrected_plus_dyntopo_SW.nc files in "
                f"{args.combined_dir}. Run `python "
                f"add_dyntopo_to_corrected_scotese.py --dyntopo-dir <...>` "
                f"first.")
        print(f"auto-discovered {len(args.ages)} ages: "
              f"{args.ages[0]}..{args.ages[-1]} Ma")
    else:
        # Validate that requested ages have combined NetCDFs
        missing = [a for a in args.ages
                   if not (args.combined_dir
                           / COMBINED_FNAME_FMT.format(age=a)).exists()]
        if missing:
            raise SystemExit(
                f"ERROR: no combined NetCDF for ages {missing} in "
                f"{args.combined_dir}. Run `python "
                f"add_dyntopo_to_corrected_scotese.py --dyntopo-dir <...>` "
                f"first.")

    make_figure(args.ages, args.combined_dir,
                figsize=tuple(args.figsize),
                write_png=not args.no_png,
                write_pdf=not args.no_pdf,
                write_csv=not args.no_csv,
                dpi=args.dpi)


if __name__ == "__main__":
    main()
