"""
=============================================================================
build_summary_stats_scotese.py  —  Per-slice & global stats (S&W workflow)
=============================================================================

PURPOSE
    Quantify how the geochem-assimilated maps differ from the original
    Scotese & Wright PaleoDEMs, **as a function of geological age**.
    The headline insight is the time evolution: when are corrections
    large, when small; which provinces drive them; how does hypsometry
    change through time.

OUTPUTS  (all written to paths_scotese.CORRECTED_DIR)

  per_slice_stats_SW.csv               primary time-dependent table.
      One row per S&W age, columns:
        t_Ma, era, continent_cells,
        n_raw, n_decluster,
        bias_before_m, bias_after_m,
        rms_before_m, rms_after_m,
        bias_reduction_m, rms_reduction_m,
        land_orig_p{5,25,50,75,95,99}_m,
        land_corr_p{5,25,50,75,95,99}_m,
        delta_mean_m, delta_rms_m, delta_p95abs_m, delta_max_m

  province_delta_summary_SW.csv        one row per (age, province):
        t_Ma, era, province, cells, n_decluster,
        delta_mean_m, delta_rms_m, delta_p95abs_m, delta_max_m,
        p50_orig_m, p50_corr_m

  before_after_summary_SW.md           human-readable global dashboard
                                       (sample-residual reductions,
                                        hypsometric envelopes, Δz stats,
                                        per-province behaviour, sample
                                        density vs era)

  era_summary_SW.md                    aggregated stats per geological
                                       era (Cenozoic, Mesozoic, Paleozoic)

INPUTS
    paths_scotese.CORRECTED_DIR / "*Ma_corrected_SW.nc"   (from assimilate_scotese.py)
    paths_scotese.CSV_PATH                                  (geochem CSV)

USAGE
    cd <project>/scripts_Scotese
    python build_summary_stats_scotese.py

DEPENDENCIES
    numpy, pandas, netCDF4  (plus assimilate_scotese for sample loading)

WORKFLOW POSITION
    Run AFTER `python assimilate_scotese.py --all`.
    full_sweep_diagnostics_scotese.py consumes per_slice_stats_SW.csv.
=============================================================================
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np, pandas as pd, netCDF4 as nc

sys.path.insert(0, str(Path(__file__).resolve().parent))
import assimilate_scotese as A
from paths_scotese import CORRECTED_DIR


# ---------------------------------------------------------------------------
PER_SLICE_CSV = CORRECTED_DIR / "per_slice_stats_SW.csv"
PROV_CSV      = CORRECTED_DIR / "province_delta_summary_SW.csv"
GLOBAL_MD     = CORRECTED_DIR / "before_after_summary_SW.md"
ERA_MD        = CORRECTED_DIR / "era_summary_SW.md"

PERCENTILES = [5, 25, 50, 75, 95, 99]
PROVS = ["Continental Arc", "Orogen", "Continental Margin", "Island Arc",
         "Extended Crust", "Basin", "Platform", "Shield", "Other"]
PROV_INDEX = {p: i for i, p in enumerate(PROVS)}


def era_of(t_Ma: float) -> str:
    """Geological era for an age in Ma — standard IUGS boundaries."""
    if t_Ma < 66:   return "Cenozoic"
    if t_Ma < 252:  return "Mesozoic"
    if t_Ma < 540:  return "Paleozoic"
    return "Proterozoic"


def percentile_dict(arr, prefix):
    out = {}
    if arr.size:
        ps = np.percentile(arr, PERCENTILES)
        for q, v in zip(PERCENTILES, ps):
            out[f"{prefix}_p{q}_m"] = float(v)
    else:
        for q in PERCENTILES:
            out[f"{prefix}_p{q}_m"] = np.nan
    return out


# ---------------------------------------------------------------------------
# Per-slice loop
# ---------------------------------------------------------------------------
def per_slice():
    rows = []
    prov_rows = []
    df_geo = A.get_geochem()
    nc_files = sorted(CORRECTED_DIR.glob("*Ma_corrected_SW.nc"),
                      key=lambda f: int(f.name.split("Ma")[0]))
    for f in nc_files:
        with nc.Dataset(f) as d:
            try:
                t = int(round(float(d.target_age_Ma)))
            except AttributeError:
                # fallback: parse age from filename
                t = int(f.name.split("Ma")[0])
            # np.asarray strips any MaskedArray wrapper netCDF4 may apply,
            # which would otherwise make np.percentile/np.partition warn
            # ("'partition' will ignore the 'mask' of the MaskedArray").
            M    = np.asarray(d.variables["M_orig"][:],         dtype=float)
            Mc   = np.asarray(d.variables["M_corrected"][:],    dtype=float)
            cont = np.asarray(d.variables["continent_mask"][:], dtype=bool)
            P    = np.asarray(d.variables["province"][:])
            lat  = np.asarray(d.variables["lat"][:])
            lon  = np.asarray(d.variables["lon"][:])

        delta = Mc - M
        dcont = delta[cont]
        land_orig = M[cont & (M > 0)]
        land_corr = Mc[cont & (Mc > 0)]

        row = dict(t_Ma=t, era=era_of(t),
                   continent_cells=int(cont.sum()))

        # sample residuals — re-prepare to match the assimilation logic exactly
        sub = A.prepare_samples(df_geo, t)
        if len(sub):
            dec = A.decluster(sub, lat, lon)
            iy = dec["iy"].astype(int).values
            ix = dec["ix"].astype(int).values
            rb = dec["z"].values - M[iy, ix]
            ra = dec["z"].values - Mc[iy, ix]
            bias_before = float(np.mean(rb)); bias_after = float(np.mean(ra))
            rms_before = float(np.sqrt(np.mean(rb**2)))
            rms_after = float(np.sqrt(np.mean(ra**2)))
            row.update(
                n_raw=int(len(sub)),
                n_decluster=int(len(dec)),
                bias_before_m=bias_before, bias_after_m=bias_after,
                rms_before_m=rms_before, rms_after_m=rms_after,
                bias_reduction_m=bias_before - bias_after,
                rms_reduction_m=rms_before - rms_after,
            )
            # province breakdown
            for prov in PROVS:
                cell_mask = (P == PROV_INDEX[prov]) & cont
                if cell_mask.sum() < 1:
                    continue
                d_in = delta[cell_mask]
                m_in = M[cell_mask]; mc_in = Mc[cell_mask]
                n_dec_p = int((dec["prov"] == prov).sum())
                prov_rows.append(dict(
                    t_Ma=t, era=era_of(t), province=prov,
                    cells=int(cell_mask.sum()),
                    n_decluster=n_dec_p,
                    delta_mean_m=float(np.mean(d_in)),
                    delta_rms_m=float(np.sqrt(np.mean(d_in**2))),
                    delta_p95abs_m=float(np.percentile(np.abs(d_in), 95)),
                    delta_max_m=float(np.max(np.abs(d_in))),
                    p50_orig_m=float(np.median(m_in)),
                    p50_corr_m=float(np.median(mc_in)),
                ))
        else:
            row.update(n_raw=0, n_decluster=0,
                       bias_before_m=np.nan, bias_after_m=np.nan,
                       rms_before_m=np.nan, rms_after_m=np.nan,
                       bias_reduction_m=np.nan, rms_reduction_m=np.nan)

        row.update(percentile_dict(land_orig, "land_orig"))
        row.update(percentile_dict(land_corr, "land_corr"))
        row.update(
            delta_mean_m=float(np.mean(dcont)) if dcont.size else np.nan,
            delta_rms_m=float(np.sqrt(np.mean(dcont**2))) if dcont.size else np.nan,
            delta_p95abs_m=float(np.percentile(np.abs(dcont), 95)) if dcont.size else np.nan,
            delta_max_m=float(np.max(np.abs(dcont))) if dcont.size else np.nan,
        )
        rows.append(row)

    pd.DataFrame(rows).to_csv(PER_SLICE_CSV, index=False)
    pd.DataFrame(prov_rows).to_csv(PROV_CSV, index=False)
    print(f"wrote {PER_SLICE_CSV}\nwrote {PROV_CSV}")
    return pd.DataFrame(rows), pd.DataFrame(prov_rows)


# ---------------------------------------------------------------------------
# Dashboard markdown
# ---------------------------------------------------------------------------
def write_global_dashboard(p_df, prov_df):
    valid = p_df.dropna(subset=["bias_before_m"])
    weights = valid["n_decluster"]
    txt = []
    txt.append("# Scotese & Wright assimilation — global before/after dashboard")
    txt.append("")
    txt.append(f"Generated from {len(p_df)} slices ({int(p_df.t_Ma.min())}..{int(p_df.t_Ma.max())} Ma).  "
               f"{len(valid)} slices had geochem samples in their temporal window.")
    txt.append("")
    txt.append("## Sample-residual statistics (declustered, sample-weighted across slices)")
    txt.append("")
    for col, label in [("bias_before_m", "Bias (input)"),
                       ("bias_after_m",  "Bias (corrected)"),
                       ("rms_before_m",  "RMS (input)"),
                       ("rms_after_m",   "RMS (corrected)")]:
        v = valid[col]
        wm = float(np.average(v, weights=weights))
        txt.append(f"- **{label}** sample-weighted mean: {wm:.0f} m   "
                   f"(median per-slice {float(v.median()):.0f} m, "
                   f"p95 {float(v.quantile(0.95)):.0f} m)")
    txt.append("")
    bias_red = valid["bias_before_m"] - valid["bias_after_m"]
    rms_red = valid["rms_before_m"] - valid["rms_after_m"]
    txt.append(f"- **Mean bias reduction per slice**: {float(bias_red.mean()):.0f} m "
               f"(median {float(bias_red.median()):.0f} m)")
    txt.append(f"- **Mean RMS reduction per slice**:  {float(rms_red.mean()):.0f} m "
               f"(median {float(rms_red.median()):.0f} m)")
    txt.append("")

    txt.append("## Hypsometric envelope (medians across slices, m)")
    txt.append("")
    txt.append("| Percentile | Input median | Corrected median |")
    txt.append("|---:|---:|---:|")
    for q in PERCENTILES:
        oi = float(p_df[f"land_orig_p{q}_m"].median())
        ci = float(p_df[f"land_corr_p{q}_m"].median())
        txt.append(f"| {q}% | {oi:.0f} | {ci:.0f} |")
    txt.append("")

    txt.append("## Δz (corrected − input) — global statistics")
    txt.append("")
    txt.append(f"- Mean Δz across all slices: {float(p_df['delta_mean_m'].mean()):.0f} m  "
               "(positive ⇒ map raised)")
    txt.append(f"- Median per-slice |Δz| RMS: {float(p_df['delta_rms_m'].median()):.0f} m")
    txt.append(f"- Median per-slice |Δz| p95: {float(p_df['delta_p95abs_m'].median()):.0f} m")
    txt.append(f"- Largest absolute Δz observed: {float(p_df['delta_max_m'].max()):.0f} m "
               "(cap = 2000 m)")
    txt.append("")

    txt.append("## Δz by tectonic province (median across all slices)")
    txt.append("")
    txt.append("| Province | n slices | median cells | median |Δz| RMS (m) | median Δz mean (m) |")
    txt.append("|---|---:|---:|---:|---:|")
    for prov in PROVS:
        sub = prov_df[prov_df["province"] == prov]
        if sub.empty: continue
        txt.append(f"| {prov} | {len(sub)} | {int(sub['cells'].median())} | "
                   f"{float(sub['delta_rms_m'].median()):.0f} | "
                   f"{float(sub['delta_mean_m'].median()):.0f} |")
    txt.append("")

    GLOBAL_MD.write_text("\n".join(txt))
    print(f"wrote {GLOBAL_MD}")


def write_era_summary(p_df, prov_df):
    txt = []
    txt.append("# Scotese & Wright assimilation — era-binned summary")
    txt.append("")
    txt.append("Geological-era aggregates of the per-slice metrics.  Useful for "
               "paper tables.  Each row is the median over slices in that era, "
               "with the IQR in parentheses.")
    txt.append("")
    for era_name in ["Cenozoic", "Mesozoic", "Paleozoic"]:
        sub = p_df[p_df.era == era_name]
        if sub.empty:
            continue
        txt.append(f"## {era_name}  ({len(sub)} slices, {int(sub.t_Ma.min())}..{int(sub.t_Ma.max())} Ma)")
        txt.append("")
        for col, label in [("n_decluster",       "Declustered samples"),
                           ("bias_before_m",     "Bias before (m)"),
                           ("bias_after_m",      "Bias after (m)"),
                           ("rms_before_m",      "RMS before (m)"),
                           ("rms_after_m",       "RMS after (m)"),
                           ("bias_reduction_m",  "Bias reduction (m)"),
                           ("rms_reduction_m",   "RMS reduction (m)"),
                           ("land_orig_p99_m",   "Land p99 input (m)"),
                           ("land_corr_p99_m",   "Land p99 corrected (m)"),
                           ("delta_rms_m",       "Δz RMS (m)"),
                           ("delta_p95abs_m",    "Δz p95 (m)")]:
            v = sub[col].dropna()
            if v.empty:
                continue
            q25, q50, q75 = v.quantile([0.25, 0.5, 0.75]).tolist()
            txt.append(f"- **{label}**: {q50:.0f} (IQR {q25:.0f} … {q75:.0f})")
        txt.append("")

        # per-province within this era
        psub = prov_df[prov_df.era == era_name]
        if not psub.empty:
            txt.append(f"### Per-province Δz in the {era_name} (median across slices)")
            txt.append("")
            txt.append("| Province | n slices | median cells | median |Δz| RMS (m) | median Δz mean (m) |")
            txt.append("|---|---:|---:|---:|---:|")
            for prov in PROVS:
                ps = psub[psub.province == prov]
                if ps.empty: continue
                txt.append(f"| {prov} | {len(ps)} | {int(ps['cells'].median())} | "
                           f"{float(ps['delta_rms_m'].median()):.0f} | "
                           f"{float(ps['delta_mean_m'].median()):.0f} |")
            txt.append("")
    ERA_MD.write_text("\n".join(txt))
    print(f"wrote {ERA_MD}")


if __name__ == "__main__":
    p_df, prov_df = per_slice()
    write_global_dashboard(p_df, prov_df)
    write_era_summary(p_df, prov_df)
