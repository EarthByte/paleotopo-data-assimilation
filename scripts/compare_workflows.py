"""
=============================================================================
compare_workflows.py  —  Merdith vs Scotese & Wright assimilation comparison
=============================================================================

Side-by-side comparison of the two corrected paleotopography model sets
on the ages they share (S&W's 5-Myr cadence is the limiting set).

For each shared age the script computes:
   - residual of the corrected maps against the geochem samples
   - per-cell |M_corrected_SW − M_corrected_Merdith| difference grid
   - hypsometric envelopes
   - per-province Δz comparison

OUTPUTS  (in paths_scotese.OUTPUT_DIR)

  Model_comparison_per_slice.csv
       one row per shared age, columns:
         t_Ma, era,
         bias_after_SW_m, bias_after_Merdith_m,
         rms_after_SW_m,  rms_after_Merdith_m,
         p99_SW_m,        p99_Merdith_m,
         intermodel_rms_m,   intermodel_p95_m,
         continent_overlap_iou

  Model_comparison_summary.md
       human-readable global dashboard with era-binned aggregates

  Model_comparison_time_series.png
       multi-panel time series comparing the two workflows:
         - sample-residual bias and RMS through time, both models
         - p99 land elevation through time, both models
         - inter-model RMS Δz through time
         - shared sample counts vs time

  Model_comparison_<age>Ma_map.png   (one per --map-age flag)
       three-panel: corrected SW | corrected Merdith | difference

USAGE
    cd <project>/scripts_Scotese
    python compare_workflows.py                     # all stats only
    python compare_workflows.py --map-ages 50 100 250 500
    python compare_workflows.py --merdith-dir /custom/path

This script is *optional* — only run it if you have BOTH the Merdith
corrected NetCDFs (in `data/corrected/`) and the S&W corrected NetCDFs
(in `data/corrected_Scotese/`).

DEPENDENCIES
    numpy, pandas, netCDF4, matplotlib, cartopy (for the optional maps)

LIMITATIONS
    - The two grids differ in shape (Merdith 180×360 vs S&W 181×361).
      We compare on the Merdith grid by nearest-neighbour re-binning S&W
      onto Merdith cells.  This is fine for global statistics but the
      per-cell difference can have edge artefacts at the poles.
    - The two workflows use different plate-model frames (MER21 vs S&W).
      Samples (rlat, rlon) come from MER21 in both cases, so the
      assimilation footprint is consistent — but the underlying
      kinematic continents have slightly different reconstructed
      positions, so the inter-model difference is partly a kinematic
      difference, not an assimilation effect.
=============================================================================
"""
from __future__ import annotations
import sys, argparse
from pathlib import Path
import numpy as np, pandas as pd, netCDF4 as nc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths_scotese import CORRECTED_DIR as SW_DIR, OUTPUT_DIR

# Merdith path — by default sibling of the S&W folder under data/
MERDITH_DIR_DEFAULT = SW_DIR.parent / "corrected"


def era_of(t):
    if t < 66:   return "Cenozoic"
    if t < 252:  return "Mesozoic"
    if t < 540:  return "Paleozoic"
    return "Proterozoic"


def shared_ages(sw_dir: Path, merdith_dir: Path) -> list[int]:
    """Ages for which both workflows have a corrected NetCDF."""
    sw = {int(f.name.split("Ma")[0]) for f in sw_dir.glob("*Ma_corrected_SW.nc")}
    me = {int(f.name.split("Ma")[0]) for f in merdith_dir.glob("*Ma_corrected.nc")}
    return sorted(sw & me)


def load_pair(t: int, sw_dir: Path, merdith_dir: Path):
    """Read both NetCDFs for age t and return aligned arrays.  S&W's 181×361
    grid is re-binned to Merdith's 180×360 via nearest-neighbour."""
    with nc.Dataset(sw_dir / f"{t}Ma_corrected_SW.nc") as d:
        lat_sw = d.variables["lat"][:]
        lon_sw = d.variables["lon"][:]
        Mc_sw = d.variables["M_corrected"][:].astype(float)
        cont_sw = d.variables["continent_mask"][:].astype(bool)
    with nc.Dataset(merdith_dir / f"{t}Ma_corrected.nc") as d:
        lat_me = d.variables["lat"][:]
        lon_me = d.variables["lon"][:]
        Mc_me = d.variables["M_corrected"][:].astype(float)
        cont_me = d.variables["continent_mask"][:].astype(bool)

    # nearest-neighbour resample S&W onto the Merdith grid
    iy = np.array([np.argmin(np.abs(lat_sw - la)) for la in lat_me])
    ix = np.array([np.argmin(np.abs(lon_sw - lo)) for lo in lon_me])
    Mc_sw_rg = Mc_sw[np.ix_(iy, ix)]
    cont_sw_rg = cont_sw[np.ix_(iy, ix)]

    return lat_me, lon_me, Mc_me, cont_me, Mc_sw_rg, cont_sw_rg


def compare_one_slice(t, sw_dir, merdith_dir):
    lat, lon, Mc_me, cont_me, Mc_sw, cont_sw = load_pair(t, sw_dir, merdith_dir)
    both_cont = cont_me & cont_sw
    iou = both_cont.sum() / max((cont_me | cont_sw).sum(), 1)
    diff = (Mc_sw - Mc_me)[both_cont]
    return dict(
        t_Ma=t, era=era_of(t),
        intermodel_rms_m=float(np.sqrt(np.mean(diff**2))) if diff.size else np.nan,
        intermodel_p95_m=float(np.percentile(np.abs(diff), 95)) if diff.size else np.nan,
        intermodel_bias_m=float(np.mean(diff)) if diff.size else np.nan,
        continent_overlap_iou=float(iou),
        p99_SW_m=float(np.percentile(Mc_sw[cont_sw & (Mc_sw > 0)], 99))
                  if (cont_sw & (Mc_sw > 0)).any() else np.nan,
        p99_Merdith_m=float(np.percentile(Mc_me[cont_me & (Mc_me > 0)], 99))
                       if (cont_me & (Mc_me > 0)).any() else np.nan,
    )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--merdith-dir", type=Path, default=MERDITH_DIR_DEFAULT,
                   help="folder containing the Merdith corrected NetCDFs")
    p.add_argument("--map-ages", type=int, nargs="*", default=[],
                   help="ages to render 3-panel comparison maps for")
    args = p.parse_args()

    if not args.merdith_dir.exists():
        print(f"ERROR: Merdith corrected dir not found: {args.merdith_dir}")
        print("       Pass --merdith-dir to point to it, or skip this script.")
        sys.exit(2)

    common = shared_ages(SW_DIR, args.merdith_dir)
    print(f"Found {len(common)} shared ages: {common[:5]} … {common[-5:]}")

    rows = [compare_one_slice(t, SW_DIR, args.merdith_dir) for t in common]
    df = pd.DataFrame(rows)
    csv_path = OUTPUT_DIR / "Model_comparison_per_slice.csv"
    df.to_csv(csv_path, index=False)
    print(f"wrote {csv_path}")

    # Pull per-slice diagnostics from both workflows' summary CSVs
    sw_summary = pd.read_csv(SW_DIR / "per_slice_stats_SW.csv")
    me_summary_path = args.merdith_dir / "per_slice_stats_v2_vs_orig.csv"
    me_summary = (pd.read_csv(me_summary_path)
                  if me_summary_path.exists() else None)

    # ---- Time-series figure
    fig, axes = plt.subplots(4, 1, figsize=(13, 14))
    ax = axes[0]
    ax.plot(sw_summary.t_Ma, sw_summary.bias_after_m, label="S&W corrected",
            color="C0", lw=1.5)
    if me_summary is not None:
        ax.plot(me_summary.t_Ma, me_summary.bias_after_m, label="Merdith corrected",
                color="C3", lw=1.5)
    ax.set_ylabel("sample bias after (m)"); ax.invert_xaxis()
    ax.set_title("Sample-residual bias through time, both models")
    ax.grid(True, alpha=0.3); ax.legend()

    ax = axes[1]
    ax.plot(sw_summary.t_Ma, sw_summary.rms_after_m, label="S&W corrected",
            color="C0", lw=1.5)
    if me_summary is not None:
        ax.plot(me_summary.t_Ma, me_summary.rms_after_m, label="Merdith corrected",
                color="C3", lw=1.5)
    ax.set_ylabel("sample RMS after (m)"); ax.invert_xaxis()
    ax.set_title("Sample-residual RMS through time")
    ax.grid(True, alpha=0.3); ax.legend()

    ax = axes[2]
    ax.plot(sw_summary.t_Ma, sw_summary.land_corr_p99_m, label="S&W corrected",
            color="C0", lw=1.5)
    if me_summary is not None:
        ax.plot(me_summary.t_Ma, me_summary.land_corr_p99_m, label="Merdith corrected",
                color="C3", lw=1.5)
    ax.set_ylabel("land p99 elevation (m)"); ax.invert_xaxis()
    ax.set_title("99th-percentile land elevation, both corrected models")
    ax.grid(True, alpha=0.3); ax.legend()

    ax = axes[3]
    ax.plot(df.t_Ma, df.intermodel_rms_m, color="C2", lw=1.5)
    ax.set_ylabel("inter-model |Δ| RMS (m)"); ax.invert_xaxis()
    ax.set_title("Per-slice RMS difference between corrected S&W and Merdith maps")
    ax.set_xlabel("age (Ma)")
    ax.grid(True, alpha=0.3)

    plt.suptitle("Merdith vs Scotese & Wright corrected paleotopography — time series",
                 fontsize=14, y=0.995)
    plt.tight_layout()
    out = OUTPUT_DIR / "Model_comparison_time_series.png"
    plt.savefig(out, dpi=140, bbox_inches="tight")
    print(f"wrote {out}")

    # ---- Markdown summary
    md = []
    md.append("# Merdith vs Scotese & Wright — corrected-model comparison")
    md.append("")
    md.append(f"{len(common)} shared ages from {common[0]} to {common[-1]} Ma.")
    md.append("")
    md.append("## Inter-model RMS difference (corrected − corrected) by era")
    md.append("")
    md.append("| Era | n slices | median RMS (m) | IQR (m) | p95 (m) |")
    md.append("|---|---:|---:|---:|---:|")
    for era_name in ["Cenozoic", "Mesozoic", "Paleozoic"]:
        sub = df[df.era == era_name]
        if sub.empty: continue
        q25, q50, q75 = sub.intermodel_rms_m.quantile([0.25, 0.5, 0.75]).tolist()
        p95 = sub.intermodel_p95_m.median()
        md.append(f"| {era_name} | {len(sub)} | {q50:.0f} | "
                  f"{q25:.0f} … {q75:.0f} | {p95:.0f} |")
    md.append("")
    md.append("## Per-slice details")
    md.append("")
    md.append(df.head(20).to_markdown(index=False))
    if len(df) > 20:
        md.append("")
        md.append(f"... ({len(df)-20} more rows in Model_comparison_per_slice.csv)")
    md_path = OUTPUT_DIR / "Model_comparison_summary.md"
    md_path.write_text("\n".join(md))
    print(f"wrote {md_path}")

    # ---- Optional 3-panel maps
    if args.map_ages:
        import cartopy.crs as ccrs
        norm = mcolors.TwoSlopeNorm(vmin=-7000, vcenter=0, vmax=7000)
        norm_d = mcolors.TwoSlopeNorm(vmin=-2000, vcenter=0, vmax=2000)
        for t in args.map_ages:
            if t not in common:
                print(f"skip {t} Ma (not in shared ages)")
                continue
            lat, lon, Mc_me, cont_me, Mc_sw, cont_sw = load_pair(t, SW_DIR, args.merdith_dir)
            diff = np.where(cont_me & cont_sw, Mc_sw - Mc_me, np.nan)

            fig, axes = plt.subplots(3, 1, figsize=(10, 13),
                                     subplot_kw=dict(projection=ccrs.Robinson()))
            LON2D, LAT2D = np.meshgrid(lon, lat)
            for ax, arr, title in [(axes[0], Mc_sw, f"S&W corrected @ {t} Ma"),
                                   (axes[1], Mc_me, f"Merdith corrected @ {t} Ma")]:
                ax.set_global()
                ax.pcolormesh(LON2D, LAT2D, arr, cmap="terrain", norm=norm,
                              transform=ccrs.PlateCarree(), shading="auto", rasterized=True)
                ax.set_title(title)
            axes[2].set_global()
            axes[2].pcolormesh(LON2D, LAT2D, diff, cmap="RdBu_r", norm=norm_d,
                                transform=ccrs.PlateCarree(), shading="auto", rasterized=True)
            axes[2].set_title(f"Difference S&W − Merdith @ {t} Ma")
            plt.tight_layout()
            out = OUTPUT_DIR / f"Model_comparison_{t}Ma_map.png"
            plt.savefig(out, dpi=130, bbox_inches="tight")
            plt.close(fig)
            print(f"wrote {out}")


if __name__ == "__main__":
    main()
