"""
=============================================================================
assimilate_scotese.py  —  Production paleo-elevation assimilation pipeline
                          for the Scotese & Wright 2018 PaleoDEMs
=============================================================================

Province-wise CDF rescaling + smoothed-residual kernel applied to the
Scotese & Wright (2018) PaleoDEM series.

WHAT THIS DOES
    For each requested age t (Ma, integer), reads the S&W PaleoDEM at the
    nearest available 5 Myr slice, finds geochem samples whose
    recon_age_Ma falls within an adaptive temporal window of t, and
    produces a corrected elevation grid that combines the S&W kinematic
    pattern with the geochem amplitude information using:

        1.  Province-wise hypsometric (CDF) rescaling within tectonic
            provinces that have ≥ N_MIN_P declustered samples.
        2.  A short-length-scale (150 km) Gaussian residual kernel that
            adjusts the rescaled field at sample-cluster locations,
            province-masked so kernels don't bleed across boundaries.

INPUTS
    - S&W PaleoDEM NetCDFs at paths_scotese.DEM_DIR  (5 Myr cadence,
      ages 0..540 Ma with two gaps: 385 and 390 Ma are missing).
    - The whole-rock igneous-rock geochemistry CSV (see data/README.md).
    - The S&W plate model files (rotation + topologies + plate polygons).

OUTPUTS  (per slice — paths_scotese.CORRECTED_DIR)
    <age>Ma_corrected_SW.nc with variables:
       lat, lon                    1-D grid axes (paleo-coordinates, S&W
                                   convention: edge-aligned −90..+90 / −180..+180)
       M_orig (lat, lon)           Scotese & Wright input elevation, m
       M_corrected (lat, lon)      assimilated elevation, m
       delta (lat, lon)            corrected − orig, m
       province (lat, lon)         int8 Tecto_Prov index (-1 = ocean)
       continent_mask (lat, lon)   bool DEM-based mask (z > COB_PROXY_M)
       n_eff (lat, lon)            effective number of samples behind cell

    A summary CSV is also incrementally written:
       paths_scotese.CORRECTED_DIR / "corrected_grids_SW_summary.csv"

CONFIGURATION  (constants defined here)
    DZ_MAX           = 2000 m     per-cell |Δz| cap
    ZMIN, ZMAX       = -9 / +7 km hard floor and ceiling (S&W goes deeper)
    RESID_LS_KM      = 150 km     Gaussian length scale
    N_MIN_P          = 30         min declustered samples per province
    KAPPA_QUALITY    = 300 m/km   σ from Missing_Risk_DeltaT
    W_MISSING_OVER40 = 0.25       down-weight low-quality samples

ADAPTIVE TEMPORAL HALF-WINDOW Δt(t) — keyed to sample sparsity, NOT to
the S&W map cadence:
    Δt(t) =  5 Myr  for 0   ≤ t < 200
    Δt(t) = 10 Myr  for 200 ≤ t < 500
    Δt(t) = 20 Myr  for 500 ≤ t < 800
    Δt(t) = 30 Myr  for 800 ≤ t ≤ 1000

USAGE
    cd <project>/scripts
    python assimilate_scotese.py 50           # single S&W slice
    python assimilate_scotese.py 50 100 200   # any list
    python assimilate_scotese.py --all        # full sweep over all S&W ages
    python assimilate_scotese.py --all --force

DEPENDENCIES
    numpy, pandas, scipy, netCDF4, pygplates

CAVEATS
    1. S&W topologies only resolve 0–100 Ma.  For older slices the
       "Continental Arc" province class is empty and the classifier
       falls back to elevation+roughness only.
    2. The geochem CSV's `rlat/rlon` columns were computed against a
       different plate model this is virtually
       indistinguishable from S&W; mismatches grow at deep time.  To
       reconstruct samples via S&W rotations, you'd need to re-map the
       PlateID column to S&W's plate IDs — see WORKFLOW_Scotese.md.

PROVENANCE / FURTHER READING
    docs/WORKFLOW_Scotese.md
    docs/Assimilation_methodology.md   (full algorithmic description)
=============================================================================
"""
from __future__ import annotations
import os, sys, json, time, argparse
import numpy as np, pandas as pd, netCDF4 as nc
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths_scotese import (CSV_PATH, CORRECTED_DIR, OUTPUT_DIR,
                           ensure_output_dirs)
from sw_io import (available_ages as sw_available_ages,
                   load_grid as sw_load_grid,
                   nearest_cell_index)
from plate_model_utils_scotese import (
    cob_mask, distance_to_subduction, subduction_zones,
    province_grid, PROV_LIST, PROV_INDEX,
    ARC_DIST_KM, MARGIN_DIST_KM, COB_PROXY_M,
)

ensure_output_dirs()
OUT_NC_DIR = CORRECTED_DIR

# ---------------------------------------------------------------------------
# Recommended config
# ---------------------------------------------------------------------------
DZ_MAX = 2000.0
ZMIN, ZMAX = -11000.0, 11000.0     # safety bounds.  Bumped ZMAX from 9000 →
                                   # 11000 m because the S&W18 PaleoDEM has
                                   # a small number of cells (~12) at Antarctic
                                   # paleo-latitudes with M_orig > 9000 m
                                   # (up to +10200 m).  The old +9000 m ZMAX
                                   # clipped them and produced 1100 m
                                   # step-function artefacts in the Δz field
                                   # against their unclipped neighbours.
RESID_LS_KM = 150.0      # per-sample residual-kernel σ (km).  Narrow
                         # kernel keeps local sample-driven structure
                         # tight; the smoothing of `r` after deposition
                         # (in assimilate_one) blends adjacent samples
                         # at a larger scale (DELTA_SMOOTH_SIGMA_FINAL).
KAPPA_QUALITY = 300.0
DELTA_SMOOTH_SIGMA       = 1.5      # first  pass on Δz: σ in cells (~165 km)
DELTA_SMOOTH_SIGMA_FINAL = 4.5      # final  pass on Δz: σ in cells (~495 km)
                                    # Moderate.  Down from 6.0 — earlier
                                    # σ=6 bled corrections too far into
                                    # neighbouring unsampled regions.
                                    # σ=4.5 still bridges the taper-
                                    # induced neutral band (~4 cells)
                                    # at province edges, but with a
                                    # tighter falloff to keep
                                    # corrections more local.
FINAL_LAND_SMOOTH_SIGMA  = 0.0      # DISABLED — the M_final-subaerial-only
                                    # smoothing was rejected (it blurred
                                    # large continent tracts uniformly
                                    # while leaving narrow artefacts in
                                    # place).  Re-enabled by setting > 0.
DELTA_MEDIAN_WINDOW = 5             # 5x5 median-filter despeckle on Δz.
                                    # Bumped from 3 → 5 to suppress speckle
                                    # CLUSTERS (up to ~12 pixels), not just
                                    # isolated single outliers.

# Cross-boundary Gaussian on Δz applied IMMEDIATELY after the
# province-edge taper, in assimilate_one() (before the residual kernel).
# Unlike the final-pass smoothing in finalise(), this one operates on
# the (M3 - M) field that's just had the taper applied, and is the main
# bridge across the thin neutral band the taper leaves at province
# boundaries.  It is masked to the continent only (not to provinces),
# so it freely crosses province boundaries.  Set to 0 to disable.
CROSS_BOUNDARY_SIGMA = 3.0          # σ in cells (~330 km).  Moderately
                                    # reduced from 4.0 to keep
                                    # corrections more local while still
                                    # bridging the post-taper neutral band
                                    # at province boundaries.

# Province-edge taper.  Each province's CDF rescaling correction is
# multiplied by a per-cell weight that fades from 1 (deep in the
# province interior) to 0 AT THE PROVINCE EDGE.  This eliminates
# the step-function discontinuity at province boundaries: both sides
# of any boundary have weight=0, so adjacent provinces' corrections
# meet at zero rather than at two different non-zero values.
#
# Weight formula:
#     w(d) = clip((d - 1) / PROV_TAPER_CELLS, 0, 1)
# where d is the cell-distance to the nearest cell of a different
# province.  Subtracting 1 from d ensures the immediate-boundary cells
# (d=1) get weight 0 exactly; cells PROV_TAPER_CELLS+1 deep get
# weight 1; smooth linear ramp in between.
#
# Set PROV_TAPER_CELLS = 0 to disable.
CONT_MASK_CLOSING_ITER = 6          # Morphological closing iterations applied
                                    # to the continent_mask in assimilate_one
                                    # to fill 1-2 cell gaps between adjacent
                                    # reconstructed polygons.  These gaps
                                    # cause the visible "white corridors"
                                    # and isolated zero-correction speckle
                                    # inside otherwise-corrected regions
                                    # (Rocky Mountains, East Asia).
                                    # Each iteration dilates then erodes by 1.
                                    # Set to 0 to disable.

PROV_TAPER_CELLS = 2                # ~220 km fade distance.  Kept SMALL
                                    # so narrow province shapes aren't
                                    # killed by the taper.  Only the very
                                    # edge cells (d=1) and the next ring
                                    # (d=2) get damped; everything d>=3
                                    # gets full correction.  The cross-
                                    # boundary smoothing (below) does the
                                    # job of bridging adjacent provinces.
W_MISSING_OVER40 = 0.25
# Floor for cells that were subaerial (M_orig > 0) in the kinematic
# prior — the corrected field must keep them above sea level so the
# assimilation never spuriously floods land that S&W18 establishes
# as subaerial.  1 m is small enough to leave near-shore lowlands
# essentially unchanged but large enough to avoid ambiguous
# zero-elevation cells in the rendered map.
SUBAERIAL_FLOOR_M = 1.0

# --- Young-age correction taper --------------------------------------------
# At 0 Ma the assimilation is bypassed entirely (M_corrected = M_orig): the
# present-day DEM is observed directly and must not be altered.  But the
# geochem correction is at full amplitude the instant one steps off 0 Ma
# (delta ~ +100 m mean over land at 5 Ma), so a chained / time-stepping
# consumer of the corrected DEMs (e.g. a stepwise landscape-evolution model)
# sees a broad continental step appear across the 0<->5 Ma boundary that is
# an artefact of the bypass, not tectonics.  The method was designed to
# amplitude-correct the kinematic prior in DEEP time, where the sparse
# geochem samples are spatially bled into a geologically plausible field;
# the youngest Cenozoic (plates near modern positions, S&W already close to
# truth) never needed it, and the abrupt 0-Ma switch was an unintended
# side effect.
#
# We therefore taper the correction smoothly to zero as t -> 0 over a short
# young-age window: ramp(t) = min(1, t / YOUNG_TAPER_T_FULL).  This leaves
# 0 Ma untouched (ramp = 0, i.e. the existing bypass), leaves ages
# >= YOUNG_TAPER_T_FULL untouched (ramp = 1, full correction), and blends
# partial corrections in between so no single step carries the full jump.
# It is a per-age, same-frame rescale of the (M_final - M_orig) field, so it
# is independent of plate motion.  Set YOUNG_TAPER_T_FULL = 0.0 to disable.
YOUNG_TAPER_T_FULL = 15.0

# Depth at which corrections fade to zero in submerged cells.  Replaces
# the previous binary "still-submerged → revert to M_orig" guard, which
# created a step-function discontinuity at the M_orig=0 coastline
# contour.  With a smooth linear fade:
#     fade(M_orig) = clip(1 + M_orig / DEPTH_FADE_M, 0, 1)
# subaerial cells get full correction (fade=1), cells at the seafloor
# down to ~DEPTH_FADE_M depth get smoothly reduced corrections, and
# cells deeper than DEPTH_FADE_M get zero correction (suppressing
# polygon-shadow drift while avoiding visible step boundaries).
DEPTH_FADE_M = 500.0

# ---------- Methodological refinements (all toggleable) ----------
# (1) Temporal kernel: weight samples by temporal distance from t.
#     "uniform"     — every sample in the window weighted equally
#     "triangular"  — w = max(0, 1 − |Age−t|/Δt(t))   ← default; smooth window edge
#     "gaussian"    — w = exp(−0.5 ((Age−t) / (Δt(t)/2))²)
TEMPORAL_KERNEL = "triangular"

# (2) Per-province "minimum" and "full-trust" sample thresholds for CDF rescaling.
#     "N" here means Kish effective sample size — n_eff = (Σw)² / Σw² —
#     applied to the declustered, temporally-weighted, optionally-pooled
#     observations.  This is robust to both per-sample uncertainty and
#     temporal-kernel down-weighting.
#
#     Behaviour as a function of n_eff:
#         n_eff <  N_MIN_P      → skip CDF rescaling for this province
#         N_MIN_P ≤ n_eff < N_FULL_P
#                               → linear shrinkage of CDF correction
#                                 (zero at N_MIN_P, full at N_FULL_P)
#         n_eff ≥ N_FULL_P      → full-amplitude CDF correction
#
#     Defaults: N_MIN_P=3 means we require at least the equivalent of three
#     independent observations to define a province's elevation
#     distribution at a given slice — fewer is too few to distinguish
#     a real distribution from an outlier and would overfit.  N_FULL_P=30
#     is the threshold above which the empirical CDF is treated as
#     fully informative (consistent with the original POC sensitivity
#     analysis on the 50 Ma slice).
N_MIN_P  = 3
N_FULL_P = 30

# (3) Province pooling: when a province has few samples, optionally pool
#     with related provinces (down-weighted) to stabilise the CDF.
PROVINCE_POOLING = True
POOL_DOWNWEIGHT  = 0.3       # weight applied to pooled non-native samples
PROVINCE_POOLS = {
    "Continental Arc":    ["Continental Arc", "Island Arc", "Orogen"],
    "Island Arc":         ["Island Arc", "Continental Arc"],
    "Continental Margin": ["Continental Margin", "Extended Crust"],
    "Basin":              ["Basin", "Continental Margin"],
    "Extended Crust":     ["Extended Crust", "Continental Margin"],
    "Orogen":             ["Orogen"],
    "Platform":           ["Platform", "Shield"],
    "Shield":             ["Shield", "Platform"],
    "Other":              ["Other"],
}

# (5) Province classifier backend.  Two options:
#       "geomorphic" — original 9-category classifier from
#                      plate_model_utils_scotese.province_grid().
#                      Continental cells are classified by elevation +
#                      roughness rules with a sample-halo overlay; sample
#                      labels come from the geochem CSV's Tecto_Prov column.
#       "hasterok"   — 15-category classifier from the
#                      Hasterok et al. (2022, ESR) global geological-
#                      province shapefile.  Polygons reconstructed per-
#                      slice via S&W rotations; samples re-classified by
#                      point-in-polygon at present-day.  Requires the
#                      shapefile at data/Hasterok_plates_provinces/shp/.
#
# Default keeps the existing scheme; set to "hasterok" to A/B compare.
PROVINCE_CLASSIFIER = "geomorphic"   # or "hasterok"

# (4) Spatial-support gate on the CDF rescaling correction.
#     The province-wise CDF rescaling builds one province-wide
#     distribution from all samples in the province globally; without
#     this gate it applies the rescaling to every cell in the province,
#     even cells > 2000 km from the nearest actual sample.  We shrink
#     the rescaling correction by a smooth distance kernel:
#
#         w(d) = clip(1 − d / LOCAL_SUPPORT_KM, 0, 1)
#
#     where d is the great-circle distance to the nearest declustered
#     sample.  Cells with d ≥ LOCAL_SUPPORT_KM receive exactly zero
#     CDF rescaling correction (the kinematic prior is preserved); the
#     smoothed-residual kernel (150 km) handles local fine structure
#     on top.  Set LOCAL_SUPPORT_KM = None or 0 to disable.
LOCAL_SUPPORT_KM = 1000.0


def dt_half(t):
    if t < 200:  return 5.0
    if t < 500:  return 10.0
    if t < 800:  return 20.0
    return 30.0


# ---------------------------------------------------------------------------
# Active province-classifier resolution.  Switch backends by setting
# PROVINCE_CLASSIFIER at the top of this file.
# ---------------------------------------------------------------------------
_CLASSIFIER_CACHE = None
def _active_classifier():
    """Return (PROV_LIST, PROV_INDEX, POOLS, grid_fn, sample_classifier).

    - grid_fn(t, lat, lon, cont, dec, M)  → 2-D int8 province grid
    - sample_classifier(lat_arr, lon_arr) → 1-D array of prov-name strings,
      or None if samples should use their geochem-CSV Tecto_Prov column.

    Result is cached after the first call.
    """
    global _CLASSIFIER_CACHE
    if _CLASSIFIER_CACHE is not None:
        return _CLASSIFIER_CACHE
    if PROVINCE_CLASSIFIER == "hasterok":
        from hasterok_provinces import (
            HASTEROK_PROV_LIST as PL,
            HASTEROK_PROV_INDEX as PI,
            HASTEROK_POOLS as POOLS_H,
            province_grid_hasterok,
            classify_samples_hasterok,
        )
        def grid_fn(t, lat, lon, cont, dec, M):
            # The Hasterok rasterisation already returns a complete
            # province raster; `dec` is unused but kept in the signature
            # for API symmetry with the geomorphic backend.
            return province_grid_hasterok(t, lat, lon, cont)
        _CLASSIFIER_CACHE = (PL, PI, POOLS_H, grid_fn, classify_samples_hasterok)
    else:
        def grid_fn(t, lat, lon, cont, dec, M):
            return province_grid(t, lat, lon, cont, dec, M)
        _CLASSIFIER_CACHE = (PROV_LIST, PROV_INDEX, PROVINCE_POOLS, grid_fn, None)
    return _CLASSIFIER_CACHE


def _R_EARTH_KM():
    return 6371.0


def nearest_sample_distance_grid(lat, lon, dec):
    """Return a (nlat, nlon) array of great-circle distances (km) from
    every grid cell to its nearest declustered sample in `dec`.

    Used by the LOCAL_SUPPORT_KM gate to suppress province-CDF
    rescaling far from any actual data.
    """
    LON, LAT = np.meshgrid(np.radians(lon), np.radians(lat))
    cos_lat = np.cos(LAT)
    sin_lat = np.sin(LAT)
    dmin = np.full(LON.shape, np.inf, dtype=np.float64)
    if dec is None or len(dec) == 0:
        return dmin
    s_lat = np.radians(dec["rlat"].values)
    s_lon = np.radians(dec["rlon"].values)
    s_sin = np.sin(s_lat); s_cos = np.cos(s_lat)
    for slat_sin, slat_cos, slon in zip(s_sin, s_cos, s_lon):
        cosD = sin_lat * slat_sin + cos_lat * slat_cos * np.cos(LON - slon)
        np.minimum(dmin, np.arccos(np.clip(cosD, -1.0, 1.0)), out=dmin)
    return _R_EARTH_KM() * dmin


# 6 geochem isostatic-elevation models in the CSV
ELEV_MODELS_MID = ["Isostatic_Elevation_absolute_km",
                   "Brown_Isostatic_Elevation_absolute_km",
                   "Davis_Isostatic_Elevation_absolute_km",
                   "Condie_Isostatic_Elevation_absolute_km",
                   "Herz_0.08_Isostatic_Elevation_absolute_km",
                   "Herz_0.38 Isostatic_Elevation_absolute_km"]
ELEV_MODELS_LOW = ["Elev_abs_low_km",
                   "Brown_Elev_abs_low_km",
                   "Davis_Elev_abs_low_km",
                   "Condie_Elev_abs_low_km",
                   "Herz_0.08_Elev_abs_low_km",
                   "Herz_0.38 Elev_abs_low_km"]
ELEV_MODELS_HIGH = ["Elev_abs_high_km",
                    "Brown_Elev_abs_high_km",
                    "Davis_Elev_abs_high_km",
                    "Condie_Elev_abs_high_km",
                    "Herz_0.08_Elev_abs_high_km",
                    "Herz_0.38 Elev_abs_high_km"]


_csv_cache = None
_PLATE_ID_CACHE = CSV_PATH.parent / "sample_plate_ids_SW.npy"

def get_geochem():
    """Load the geochem CSV, ensure each sample has an S&W plate ID
    assigned (cached on disk after the first run since the point-in-
    polygon test is the slowest step in the whole pipeline)."""
    global _csv_cache
    if _csv_cache is not None:
        return _csv_cache.copy()
    df = pd.read_csv(CSV_PATH, low_memory=False)

    # Try to load cached plate IDs from previous run.  Cache validity is
    # judged by row count alone (the geochem CSV is treated as immutable).
    if _PLATE_ID_CACHE.exists():
        cached = np.load(_PLATE_ID_CACHE)
        if len(cached) == len(df):
            df["sw_plate_id"] = cached
            print(f"  loaded cached S&W plate IDs from {_PLATE_ID_CACHE.name} "
                  f"({(cached != 0).sum()} continental, {(cached == 0).sum()} ocean)")
        else:
            print(f"  cached plate-ID file size mismatch — will rebuild")

    if "sw_plate_id" not in df.columns:
        print(f"  assigning S&W plate IDs by point-in-polygon "
              f"({len(df)} samples — this is a one-time ~30 s job) …")
        sr = get_reconstructor()
        pids = sr.assign_plate_ids(df)
        df["sw_plate_id"] = pids
        try:
            np.save(_PLATE_ID_CACHE, pids)
            print(f"  cached plate IDs to {_PLATE_ID_CACHE}")
        except Exception as e:
            print(f"  could not write plate-ID cache: {e}")

    # Drop oceanic samples — anchoring them at present-day would put them
    # in the ocean at older times when continents have moved away.
    n_total = len(df)
    df = df[df["sw_plate_id"] != 0].reset_index(drop=True)
    n_dropped = n_total - len(df)
    if n_dropped:
        print(f"  dropped {n_dropped} samples not contained in any "
              f"continental polygon ({100*n_dropped/n_total:.1f} %)")

    _csv_cache = df
    return _csv_cache.copy()


# ----------------------------------------------------------------------
# Per-slice sample reconstruction (correctness fix) — see WORKFLOW.
# ----------------------------------------------------------------------
from sample_reconstruct_scotese import ScoteseSampleReconstructor
_reconstructor = None
def get_reconstructor():
    global _reconstructor
    if _reconstructor is None:
        _reconstructor = ScoteseSampleReconstructor()
    return _reconstructor


# ---------------------------------------------------------------------------
# Sample preparation, declustering
# (numerically identical; tied here for self-containment)
# ---------------------------------------------------------------------------
def prepare_samples(df, t):
    """Filter samples to the ±Δt(t) window around target age t and
    RECONSTRUCT each sample's paleo-coordinates TO age t using S&W
    rotations and the plate ID assigned in get_geochem()."""
    Δ = dt_half(t)
    mask = (df["Age_Ma"] >= t - Δ) & (df["Age_Ma"] <= t + Δ)
    df = df[mask].copy()
    df = df.dropna(subset=["Lat", "Lon"])
    df = df[df.get("Missing_Risk_DeltaT (km)", 0).fillna(0) <= 5].copy()
    if df.empty:
        return df

    # --- per-slice reconstruction via S&W rotations + assigned plate IDs
    df = df.reset_index(drop=True)
    if "sw_plate_id" not in df.columns:
        # caller used get_geochem() incorrectly — assign on the fly
        df["sw_plate_id"] = get_reconstructor().assign_plate_ids(df)
    rlat_t, rlon_t = get_reconstructor().reconstruct(df, t)
    df["rlat"] = rlat_t
    df["rlon"] = rlon_t
    df = df.dropna(subset=["rlat", "rlon"]).reset_index(drop=True)
    if df.empty:
        return df
    # ------------------------------------------------------------------

    mid = df[ELEV_MODELS_MID].astype(float)
    z_mean_km = mid.mean(axis=1, skipna=True)
    z_std_km = mid.std(axis=1, skipna=True).fillna(0.0)
    low = df[ELEV_MODELS_LOW].astype(float).mean(axis=1, skipna=True)
    high = df[ELEV_MODELS_HIGH].astype(float).mean(axis=1, skipna=True)
    sigma_thick_km = ((high - low) / 4.0).abs().fillna(0.5)
    miss_risk = df.get("Missing_Risk_DeltaT (km)", pd.Series(0)).fillna(0.0)
    sigma_total_m = np.sqrt((z_std_km*1000.0)**2 +
                            (sigma_thick_km*1000.0)**2 +
                            (KAPPA_QUALITY*miss_risk)**2)
    sigma_total_m = np.clip(sigma_total_m, 200.0, 4000.0)
    df["z_obs_m"] = z_mean_km * 1000.0
    df["sigma_m"] = sigma_total_m
    w = 1.0 / (df["sigma_m"]**2)
    if "Missing_Over40pct" in df.columns:
        flag = df["Missing_Over40pct"].astype(str).str.upper().eq("TRUE")
        w = w.where(~flag, w * W_MISSING_OVER40)

    # --- Temporal-distance weight (refinement (1)) ---
    if TEMPORAL_KERNEL != "uniform":
        abs_dt = (df["Age_Ma"] - float(t)).abs()
        if TEMPORAL_KERNEL == "triangular":
            tw = np.clip(1.0 - abs_dt / float(Δ), 0.0, 1.0)
        elif TEMPORAL_KERNEL == "gaussian":
            sigma_t = max(float(Δ) / 2.0, 1e-9)
            tw = np.exp(-0.5 * (abs_dt / sigma_t) ** 2)
        else:
            raise ValueError(f"Unknown TEMPORAL_KERNEL: {TEMPORAL_KERNEL!r}")
        w = w * tw

    df["w_raw"] = w

    # Province assignment: if a non-geomorphic classifier is in use,
    # override the CSV's Tecto_Prov with the alternative scheme's
    # per-sample point-in-polygon classification (cached).
    _PL, _PI, _POOLS, _GRIDFN, _CLASSIFIER = _active_classifier()
    if _CLASSIFIER is not None:
        df["Tecto_Prov"] = _CLASSIFIER(df["Lat"].values, df["Lon"].values)
    df["Tecto_Prov"] = df["Tecto_Prov"].fillna("Other")
    return df.dropna(subset=["z_obs_m"])


def decluster(df, lat, lon):
    """Bin by 1° (lat, lon) × Tecto_Prov; replace each group with its
    weighted-median and an effective σ.  S&W grid is edge-aligned so the
    cell mapping uses nearest_cell_index from sw_io."""
    iy = nearest_cell_index(lat, df["rlat"].values)
    ix = nearest_cell_index(lon, df["rlon"].values)
    df = df.assign(iy=iy, ix=ix)
    rows = []
    for (iy_, ix_, prov), g in df.groupby(["iy", "ix", "Tecto_Prov"]):
        w = g["w_raw"].values; z = g["z_obs_m"].values; sig = g["sigma_m"].values
        order = np.argsort(z); zs = z[order]; ws = w[order]
        # Guard against an all-zero-weight bin (degenerate sample group):
        # fall back to unweighted median rather than dividing by zero.
        w_tot = float(np.sum(ws))
        if w_tot <= 0.0 or not np.isfinite(w_tot):
            z_med = float(np.median(zs))
        else:
            cum = np.cumsum(ws) / w_tot
            z_med = float(zs[np.searchsorted(cum, 0.5)])
        kish = (w.sum())**2 / (np.square(w).sum() + 1e-30)
        z_std = float(np.sqrt(max(np.var(z), 0.0)))
        sigma_eff = max(z_std, sig.mean() / max(np.sqrt(kish), 1.0))
        w_eff = w.sum() / np.sqrt(max(kish, 1.0))
        rows.append(dict(iy=int(iy_), ix=int(ix_), prov=prov,
                         rlat=float(lat[iy_]), rlon=float(lon[ix_]),
                         z=z_med, sigma=sigma_eff, w=w_eff,
                         n=len(g), n_eff=kish))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Province-wise CDF rescaling
#
# Refinements:
#   (2) per-province N thresholds: skip CDF if n_eff < N_MIN_P
#   (3) shrinkage: between N_MIN_P and N_FULL_P, the CDF correction is
#       linearly scaled from 0 to 1.  Above N_FULL_P, full correction.
#   (4) province pooling: when the native province has few samples, pool
#       with related provinces (PROVINCE_POOLS) using POOL_DOWNWEIGHT on
#       the non-native contributions.
# ---------------------------------------------------------------------------
def _samples_for_province(dec, prov, use_pooling: bool = PROVINCE_POOLING):
    """Return a DataFrame of declustered samples to use when building the
    target CDF for province `prov`.  Native samples carry their original
    weight; pooled (non-native) samples are down-weighted by POOL_DOWNWEIGHT.

    The pooling table is the active classifier's (geomorphic or
    hasterok), not the module-level PROVINCE_POOLS constant.
    """
    native = dec[dec["prov"] == prov].copy()
    if native.empty or not use_pooling:
        return native
    _PL, _PI, _POOLS, _GRIDFN, _ = _active_classifier()
    pool_names = _POOLS.get(prov, [prov])
    others = pool_names[1:] if pool_names[0] == prov else [p for p in pool_names if p != prov]
    others_df = dec[dec["prov"].isin(others)].copy()
    if others_df.empty:
        return native
    others_df["w"] = others_df["w"] * POOL_DOWNWEIGHT
    return pd.concat([native, others_df], ignore_index=True)


def _kish_n_eff(w: np.ndarray) -> float:
    w = np.asarray(w, dtype=float)
    s = w.sum()
    if s <= 0 or len(w) == 0:
        return 0.0
    return float(s * s / max((w * w).sum(), 1e-30))


def province_rescale(M, P_grid, dec, cont_mask,
                     n_min_p: int = N_MIN_P, n_full_p: int = N_FULL_P,
                     use_pooling: bool = PROVINCE_POOLING):
    """Province-wise CDF rescaling.  Iterates over the ACTIVE
    classifier's PROV_LIST (geomorphic 9-category or Hasterok 16-category)
    so the function transparently handles either backend."""
    _PL, _PI, _POOLS, _GRIDFN, _ = _active_classifier()
    Mp = M.copy()
    for prov_idx, prov in enumerate(_PL):
        cell_mask = (P_grid == prov_idx) & cont_mask
        if cell_mask.sum() < 50:
            continue

        sub = _samples_for_province(dec, prov, use_pooling=use_pooling)
        if sub.empty:
            continue

        # Effective sample count via Kish on the (possibly pooled) weights
        n_eff = _kish_n_eff(sub["w"].values)
        if n_eff < n_min_p:
            continue                            # too sparse — skip CDF for this province

        # Shrinkage factor: 0 at N_MIN_P, 1 at N_FULL_P (refinement 3)
        shrink = np.clip((n_eff - n_min_p) / max(n_full_p - n_min_p, 1.0), 0.0, 1.0)

        z_obs = sub["z"].values
        w_obs = sub["w"].values
        order = np.argsort(z_obs)
        z_sorted = z_obs[order]; w_sorted = w_obs[order]
        c_obs = np.cumsum(w_sorted); c_obs /= c_obs[-1]

        z_grid_sorted = np.sort(M[cell_mask])
        c_grid = np.linspace(0, 1, len(z_grid_sorted))

        v = M[cell_mask]
        q = np.interp(v, z_grid_sorted, c_grid)
        taper = np.minimum(np.clip(q / 0.02, 0, 1), np.clip((1 - q) / 0.02, 0, 1))
        v_remap = np.interp(q, c_obs, z_sorted)
        # Combined shrinkage + edge taper:
        new_v = v + shrink * taper * (v_remap - v)
        delta = np.clip(new_v - v, -DZ_MAX, DZ_MAX)
        Mp[cell_mask] = v + delta
    return Mp


# ---------------------------------------------------------------------------
# Localised smoothed-residual kernel (province-masked)
# ---------------------------------------------------------------------------
def smoothed_residual(M_base, P_grid, dec, lat, lon, cont_mask, ls_km=RESID_LS_KM):
    """For each declustered sample, deposit a Gaussian kernel on the grid
    around the sample's reconstructed paleo-position.  Continent-masked
    only (no province mask — that produces rectangular cell-aligned
    cut-outs around samples sitting on province boundaries).  The
    per-province CDF rescaling already handles province-distinct bulk
    corrections."""
    nlat, nlon = M_base.shape
    if dec.empty:
        return np.zeros_like(M_base), np.zeros_like(M_base)

    sample_lat = dec["rlat"].values; sample_lon = dec["rlon"].values
    iy = dec["iy"].astype(int).values; ix = dec["ix"].astype(int).values
    sample_resid = dec["z"].values - M_base[iy, ix]
    sample_w = dec["w"].values

    num = np.zeros_like(M_base); den = np.zeros_like(M_base); den_sq = np.zeros_like(M_base)
    # Bbox cutoff at 5·ℓ so the Gaussian (~ exp(-12.5)) fades to ~zero at
    # the bbox edge — no rectangular bbox artefacts.
    KERNEL_SIGMAS = 5.0
    bbox_lat_cells = int(np.ceil(KERNEL_SIGMAS * ls_km / 111.32)) + 1
    R = 6371.0

    for k in range(len(dec)):
        ic, jc = iy[k], ix[k]
        i0, i1 = max(0, ic - bbox_lat_cells), min(nlat, ic + bbox_lat_cells + 1)
        cos_lat = max(np.cos(np.radians(sample_lat[k])), 0.05)
        bbox_lon_cells = int(np.ceil(KERNEL_SIGMAS * ls_km / (111.32 * cos_lat))) + 1
        j_idx = (np.arange(jc - bbox_lon_cells, jc + bbox_lon_cells + 1)) % nlon
        lat_sub = lat[i0:i1]
        lon_sub = lon[j_idx]
        LAT2D, LON2D = np.meshgrid(lat_sub, lon_sub, indexing="ij")
        lat1 = np.radians(LAT2D); lat2 = np.radians(sample_lat[k])
        dlat = lat2 - lat1
        dlon = np.radians(sample_lon[k] - LON2D)
        a = np.sin(dlat/2)**2 + np.cos(lat1)*np.cos(lat2)*np.sin(dlon/2)**2
        d = 2 * R * np.arcsin(np.sqrt(np.clip(a, 0, 1)))
        K = np.exp(-0.5 * (d / ls_km)**2)
        K[d > KERNEL_SIGMAS * ls_km] = 0.0
        # Mask to continental cells only — no province mask.
        cont_sub = cont_mask[i0:i1, :][:, j_idx]
        K = np.where(cont_sub, K, 0.0)
        wK = K * sample_w[k]
        np.add.at(num, (slice(i0, i1), j_idx), wK * sample_resid[k])
        np.add.at(den, (slice(i0, i1), j_idx), wK)
        np.add.at(den_sq, (slice(i0, i1), j_idx), wK**2)

    valid = den > 1e-30
    r = np.where(valid, num / np.where(valid, den, 1), 0.0)
    n_eff = np.where(valid, (den**2) / (den_sq + 1e-30), 0.0)
    alpha = 1 - np.exp(-n_eff)
    r *= alpha
    return np.where(cont_mask, r, 0.0), n_eff


def province_edge_taper(P_grid, cont_mask, taper_cells: float = PROV_TAPER_CELLS):
    """Compute a per-cell weight in [0, 1] that fades each province's
    correction smoothly to ZERO at its boundary.

    For every continental cell, the weight is
        w(d) = clip((d - 1) / taper_cells, 0, 1)
    where d is the Euclidean distance (in grid cells) from this cell
    to the nearest cell of a DIFFERENT province (or outside cont_mask).

    - d = 1 (boundary cell, immediately adjacent to a different province):
      w = 0  →  no correction applied to this cell
    - d = taper_cells + 1 (deeper interior):
      w = 1  →  full correction
    - Linear ramp between

    Because cells RIGHT AT the boundary have w=0 on both sides, adjacent
    provinces' corrections meet at zero, not at two different non-zero
    values.  This eliminates the step-function discontinuity at province
    boundaries by construction.
    """
    if taper_cells is None or taper_cells <= 0:
        return np.ones_like(P_grid, dtype=np.float32)
    from scipy.ndimage import distance_transform_edt
    nlat, nlon = P_grid.shape
    weight = np.zeros((nlat, nlon), dtype=np.float32)
    unique_provs = np.unique(P_grid[cont_mask])
    for p in unique_provs:
        if int(p) < 0:
            continue
        same = (P_grid == p) & cont_mask
        if not same.any():
            continue
        d = distance_transform_edt(same)
        # Subtract 1 so cells right at the boundary (d=1) get weight 0.
        weight[same] = np.clip((d[same] - 1.0) / float(taper_cells),
                               0.0, 1.0)
    return weight


def despeckle_province_raster(P_grid, cont_mask, min_neighbours: int = 3,
                              window: int = 3, max_passes: int = 3):
    """Replace single-pixel and 2-cell-wide province artefacts with the
    local majority province.

    For each continental cell, count how many cells in the `window`×`window`
    neighbourhood share its province.  If fewer than `min_neighbours`
    (including the cell itself), reassign the cell to the most common
    province in its neighbourhood.  Iterate up to `max_passes` times so
    that 2-cell-wide bands also get absorbed.

    This eliminates the thin 1-pixel-wide stripes that appear in the
    rasterised province grid (e.g. Late Cretaceous Cordilleran accretionary
    terranes, which Hasterok represents as very narrow elongated polygons)
    and the isolated single-pixel speckle in mountainous regions.
    """
    out = P_grid.copy()
    nlat, nlon = out.shape
    pad = window // 2
    for _ in range(max_passes):
        changed = 0
        # Look up unique province codes once
        unique_provs = np.unique(out[cont_mask])
        unique_provs = unique_provs[unique_provs >= 0]
        # Build per-province occurrence count via 3x3 box filter on
        # binary indicator arrays — much faster than per-cell loop
        from scipy.ndimage import uniform_filter
        counts = {}
        for p in unique_provs:
            indicator = (out == p).astype(np.float32)
            box = uniform_filter(indicator, size=window, mode="nearest")
            counts[int(p)] = (box * (window * window)).astype(np.int16)
        # For each cell, find dominant province in its neighbourhood
        # and check whether the current province is well-supported.
        # Stack counts into (nprov, nlat, nlon)
        stack = np.stack([counts[int(p)] for p in unique_provs], axis=0)
        argmax_idx = stack.argmax(axis=0)
        majority = unique_provs[argmax_idx]
        majority_count = stack.max(axis=0)
        # Count of the cell's own province in its window
        own_count = np.zeros_like(out, dtype=np.int16)
        for p in unique_provs:
            same = (out == p)
            own_count[same] = counts[int(p)][same]
        reassign = cont_mask & (own_count < min_neighbours) & (majority != out) & (majority_count >= min_neighbours)
        out = np.where(reassign, majority, out).astype(np.int8)
        changed = int(reassign.sum())
        if changed == 0:
            break
    return out


def _masked_gaussian_smooth(arr, mask, sigma_cells):
    """Gaussian-smooth `arr` over the True region of `mask` only (no
    bleed into the False region).  Standard convolution-renormalisation
    trick: smooth (arr*mask) and (mask) separately, then divide."""
    from scipy.ndimage import gaussian_filter
    arr_in = np.where(mask, arr.astype(float), 0.0)
    smoothed_arr = gaussian_filter(arr_in, sigma=sigma_cells, mode="nearest")
    smoothed_mask = gaussian_filter(mask.astype(float), sigma=sigma_cells, mode="nearest")
    out = np.where(smoothed_mask > 1e-3, smoothed_arr / np.maximum(smoothed_mask, 1e-9), 0.0)
    return np.where(mask, out, 0.0)


def finalise(M_orig, M_cand, cont_mask):
    """Compose the final corrected elevation field with smooth boundaries.

       1. raw Δz = M_cand − M_orig  (from province CDF + residual kernel)
       2. depth-weighted fade — suppresses corrections smoothly in
          submerged cells (replaces a previous hard "revert to M_orig"
          guard that created step-function boundaries at the coastline)
       3. continent-masked Gaussian smoothing of Δz
       4. per-cell |Δz| cap at DZ_MAX
       5. global safety floor/ceiling [ZMIN, ZMAX]
       6. zero correction outside the continent mask
       7. subaerial-floor guard (M_orig > 0 ⇒ M_final ≥ SUBAERIAL_FLOOR_M)
       8. final continent-masked Gaussian smoothing pass so the cap +
          floor + fade boundaries don't show up as visible step lines.
    """
    delta = M_cand - M_orig

    # (2) Depth-weighted fade: subaerial cells get full correction;
    # cells deeper than DEPTH_FADE_M get zero correction; linear ramp
    # in between.  This suppresses polygon-shadow drift in deep ocean
    # cells while keeping the coastline transition smooth.
    if DEPTH_FADE_M and DEPTH_FADE_M > 0:
        fade = np.clip(1.0 + M_orig / DEPTH_FADE_M, 0.0, 1.0)
        delta = delta * fade

    # (2b) Median-filter despeckle of Δz inside the continent mask.
    # Replaces each cell's Δz with the median of its NxN neighbourhood,
    # killing single-pixel salt-and-pepper noise from per-sample
    # residual deposition before the Gaussian smoothing.
    if DELTA_MEDIAN_WINDOW and DELTA_MEDIAN_WINDOW >= 3:
        from scipy.ndimage import median_filter
        delta_in = np.where(cont_mask, delta, np.nan)
        # nan-safe median: convert NaN to 0 with a renormalised
        # contribution mask, or just use median_filter and reinstate
        # NaN.  Cleanest: median_filter with the cell's NaN handled
        # by setting nan→0 before; result is reasonable inside the mask.
        delta_in = np.where(cont_mask, delta, 0.0)
        delta_med = median_filter(delta_in, size=DELTA_MEDIAN_WINDOW,
                                  mode="nearest")
        delta = np.where(cont_mask, delta_med, delta)

    # (3) First continent-masked Gaussian smoothing of Δz
    if DELTA_SMOOTH_SIGMA > 0:
        delta = _masked_gaussian_smooth(delta, cont_mask, DELTA_SMOOTH_SIGMA)

    # (4) Per-cell cap, (5) global bounds, (6) zero outside continent
    delta = np.clip(delta, -DZ_MAX, DZ_MAX)
    M_final = np.clip(M_orig + delta, ZMIN, ZMAX)
    M_final = np.where(cont_mask, M_final, M_orig)

    # (7) Subaerial-floor guard
    was_subaerial = (M_orig > 0)
    M_final = np.where(was_subaerial & (M_final < SUBAERIAL_FLOOR_M),
                       SUBAERIAL_FLOOR_M, M_final)

    # (8) Optional Δz final smoothing pass (legacy; disabled by default
    # since FINAL_LAND_SMOOTH_SIGMA replaces it).
    if DELTA_SMOOTH_SIGMA_FINAL > 0:
        delta_final = M_final - M_orig
        delta_final = _masked_gaussian_smooth(delta_final, cont_mask,
                                              DELTA_SMOOTH_SIGMA_FINAL)
        M_final = M_orig + delta_final
        M_final = np.clip(M_final, ZMIN, ZMAX)
        M_final = np.where(cont_mask, M_final, M_orig)
        M_final = np.where(was_subaerial & (M_final < SUBAERIAL_FLOOR_M),
                           SUBAERIAL_FLOOR_M, M_final)

    # (9) Final smoothing pass on M_final ITSELF (the merged corrected
    # elevation), restricted to subaerial cells (M_final > 0).
    # Bathymetry (M_final ≤ 0) is preserved EXACTLY — only land gets
    # smoothed.  This addresses the visual sharp edges that come from
    # any source (Δz patches, residual-kernel speckle, OR inherited
    # M_orig features), all in one pass, while leaving bathymetric
    # detail intact.
    if FINAL_LAND_SMOOTH_SIGMA and FINAL_LAND_SMOOTH_SIGMA > 0:
        subaerial = (M_final > 0) & cont_mask
        if subaerial.any():
            M_land_smooth = _masked_gaussian_smooth(M_final, subaerial,
                                                   FINAL_LAND_SMOOTH_SIGMA)
            M_final = np.where(subaerial, M_land_smooth, M_final)
            # Smoothing may have nudged some subaerial cells just below
            # the floor — re-enforce.
            M_final = np.where(was_subaerial & (M_final < SUBAERIAL_FLOOR_M),
                               SUBAERIAL_FLOOR_M, M_final)

    return M_final


# ---------------------------------------------------------------------------
# Young-age correction taper
# ---------------------------------------------------------------------------
def young_taper_ramp(t):
    """Linear taper weight for the geochem correction near present day.

    0 at 0 Ma (the present-day bypass), rising linearly to 1 at
    YOUNG_TAPER_T_FULL, and clamped to 1 for older ages.  See the
    YOUNG_TAPER_T_FULL comment block for the rationale.
    """
    if YOUNG_TAPER_T_FULL is None or YOUNG_TAPER_T_FULL <= 0:
        return 1.0
    return float(min(1.0, max(0.0, t / YOUNG_TAPER_T_FULL)))


def apply_young_taper(t, M_orig, M_final, cont_mask):
    """Scale the assimilation correction (M_final - M_orig) by the young-age
    ramp, re-applying finalise()'s invariants (safety bounds, ocean =
    M_orig, subaerial floor).  A no-op for ages >= YOUNG_TAPER_T_FULL.

    This is a same-age, per-cell rescale of the correction field, so it
    carries no plate-motion dependence: it only changes HOW MUCH of the
    already-computed correction is applied at young ages, never WHERE.
    """
    r = young_taper_ramp(t)
    if r >= 1.0:
        return M_final
    delta = (M_final - M_orig) * r            # shrink toward the S&W prior
    Mt = np.clip(M_orig + delta, ZMIN, ZMAX)
    Mt = np.where(cont_mask, Mt, M_orig)      # ocean untouched
    was_subaerial = (M_orig > 0)
    Mt = np.where(was_subaerial & (Mt < SUBAERIAL_FLOOR_M),
                  SUBAERIAL_FLOOR_M, Mt)
    return Mt


# ---------------------------------------------------------------------------
# Per-slice orchestrator
# ---------------------------------------------------------------------------
def assimilate_one(t, save_nc=True):
    t0 = time.time()
    M, lat, lon = sw_load_grid(t)
    cont = cob_mask(t, lat, lon, M)
    # Fill 1-2 cell gaps between adjacent reconstructed polygons by
    # morphological closing.  These gaps appear as "white corridors"
    # and isolated unprocessed speckle inside otherwise-corrected
    # regions because cells with cont=False get M_final = M_orig
    # (no correction).
    if CONT_MASK_CLOSING_ITER and CONT_MASK_CLOSING_ITER > 0:
        from scipy.ndimage import binary_closing
        n_before = int(cont.sum())
        cont = binary_closing(cont, iterations=CONT_MASK_CLOSING_ITER)
        n_after = int(cont.sum())
        if n_after > n_before:
            print(f"             continent_mask: filled {n_after - n_before} "
                  f"polygon-gap cells via {CONT_MASK_CLOSING_ITER}-iter closing")

    # ---- Special case: present-day ----
    # At t=0 Ma the kinematic prior IS the modern observed topography.
    # We have direct knowledge of present-day elevations, so applying
    # any geochemically-derived correction here is both unnecessary and
    # would degrade the field.  Return M_corrected = M_orig with zero
    # delta, write the NetCDF with the standard schema, and short-circuit
    # all the per-slice diagnostics.
    if int(round(float(t))) == 0:
        print(f"[    0 Ma] bypassing assimilation — corrected = original "
              "(modern topography known directly)")
        M_final = M.copy()
        _PL_b, _PI_b, _, _, _ = _active_classifier()
        P_g = np.where(cont, _PI_b.get("Other", 0), -1).astype(np.int8)
        n_eff_grid = np.zeros_like(M)
        dec = pd.DataFrame()
        elapsed = time.time() - t0
        summary = dict(
            t_Ma=0, wall_s=elapsed,
            n_samples=0, n_decluster=0,
            continent_cells=int(cont.sum()),
            bias_before_m=0.0, bias_after_m=0.0,
            rms_before_m=0.0,  rms_after_m=0.0,
            p99_before_m=float(np.nanpercentile(M[cont & (M > 0)], 99))
                         if cont.any() else np.nan,
            p99_after_m =float(np.nanpercentile(M[cont & (M > 0)], 99))
                         if cont.any() else np.nan,
            delta_rms_m=0.0,
        )
        if save_nc:
            out = OUT_NC_DIR / f"0Ma_corrected_SW.nc"
            if out.exists():
                out.unlink()
            with nc.Dataset(out, "w") as f:
                f.createDimension("lat", len(lat))
                f.createDimension("lon", len(lon))
                f.createVariable("lat", "f4", ("lat",))[:] = lat
                f.createVariable("lon", "f4", ("lon",))[:] = lon
                f.createVariable("M_orig",      "f4", ("lat","lon"))[:] = M
                f.createVariable("M_corrected", "f4", ("lat","lon"))[:] = M_final
                f.createVariable("delta",       "f4", ("lat","lon"))[:] = M_final - M
                f.createVariable("province",    "i1", ("lat","lon"))[:] = P_g
                f.createVariable("continent_mask", "i1", ("lat","lon"))[:] = cont.astype(np.int8)
                f.createVariable("n_eff",       "f4", ("lat","lon"))[:] = n_eff_grid
                f.target_age_Ma = 0.0
                f.dataset = "Scotese & Wright 2018 PaleoDEMs"
                f.config  = "present-day bypass (M_corrected = M_orig)"
        return summary, M, M_final, P_g, cont, dec

    df = prepare_samples(get_geochem(), t)
    print(f"[{t:5.0f} Ma] samples in window ±{dt_half(t):.0f} Myr: {len(df)}")
    if not df.empty:
        # Validation: samples whose reconstructed position lands offshore
        # are a sign of plate-ID misassignment.  Drop and warn.
        #
        # The criterion is the rasterised continental polygon mask alone —
        # we deliberately do NOT additionally require M > 0 m here.  The
        # S&W PaleoDEM contains many continental cells just below sea level
        # (passive-margin shelves, rift basins) that still legitimately
        # record continental-crust geochemistry and should constrain the
        # assimilation.  The stricter `cont & (M > 0) & (Mc > 0)` criterion
        # is applied only at the visualisation layer (see render_videos_*
        # and make_pipeline_illustration_figure.py) to keep overlay points
        # off the rendered ocean.
        iy = nearest_cell_index(lat, df["rlat"].values)
        ix = nearest_cell_index(lon, df["rlon"].values)
        on_cont = cont[iy, ix]
        n_total = len(df)
        df = df[on_cont].reset_index(drop=True)
        n_oceanic = n_total - len(df)
        if n_oceanic:
            pct = 100 * n_oceanic / n_total
            warn = "  ⚠" if pct > 25 else ""
            print(f"             {n_oceanic}/{n_total} samples landed offshore after "
                  f"reconstruction → dropped ({pct:.1f} %){warn}")

    _PL, _PI, _POOLS, _GRIDFN, _ = _active_classifier()

    if df.empty:
        M_final = M.copy()
        P_g = np.where(cont, _PI.get("Other", 0), -1).astype(np.int8)
        n_eff_grid = np.zeros_like(M); dec = pd.DataFrame()
    else:
        dec = decluster(df, lat, lon)
        print(f"             declustered: {len(dec)} ({dec['prov'].value_counts().to_dict()})")
        print(f"             classifier:  {PROVINCE_CLASSIFIER}  "
              f"({len(_PL)} provinces)")
        P_g = _GRIDFN(t, lat, lon, cont, dec, M)
        # Remove 1-pixel-wide province artefacts (thin rasterised polygons,
        # isolated single-cell labels) before the CDF rescaling so the
        # rescaled field doesn't carry them through as bright stripes.
        P_g = despeckle_province_raster(P_g, cont)
        M3_raw = province_rescale(M, P_g, dec, cont, n_min_p=N_MIN_P)

        # --- Spatial-support gate on the CDF rescaling ---
        # The CDF rescaling builds one province-wide distribution from
        # all samples in the province globally; without gating it
        # applies that rescaling to every cell in the province, even
        # cells > 2000 km from the nearest actual sample (e.g.
        # Antarctica, African passive margins, the Tethys / Mediterranean
        # interior).  We shrink the rescaling correction by a smooth
        # distance kernel that reaches exactly zero at LOCAL_SUPPORT_KM.
        # The 150-km residual kernel handles fine local structure on
        # top of the gated rescaling.
        if LOCAL_SUPPORT_KM and LOCAL_SUPPORT_KM > 0:
            dist_km = nearest_sample_distance_grid(lat, lon, dec)
            w_local = np.clip(1.0 - dist_km / LOCAL_SUPPORT_KM, 0.0, 1.0)
            dz_rescale = M3_raw - M
            M3 = M + w_local * dz_rescale
            far_frac = float((dist_km > LOCAL_SUPPORT_KM)[cont].mean())
            print(f"             local-support gate: L={LOCAL_SUPPORT_KM:.0f} km; "
                  f"{100*far_frac:.0f} % of continental cells beyond "
                  f"support → zero rescaling")
        else:
            M3 = M3_raw

        # Province-edge taper: each province's rescaling correction
        # fades to zero at its boundary, so adjacent provinces' Δz
        # meet at zero rather than two different non-zero values.
        if PROV_TAPER_CELLS and PROV_TAPER_CELLS > 0:
            w_edge = province_edge_taper(P_g, cont, PROV_TAPER_CELLS)
            M3 = M + w_edge * (M3 - M)

        # Cross-boundary Gaussian smoothing of (M3 - M) — bridges the
        # thin neutral band the taper leaves at province boundaries by
        # blending corrections from adjacent provinces.  Continent-
        # masked but NOT province-masked, so it freely crosses province
        # boundaries.
        if CROSS_BOUNDARY_SIGMA and CROSS_BOUNDARY_SIGMA > 0:
            dz_prov = M3 - M
            dz_prov = _masked_gaussian_smooth(dz_prov, cont,
                                              CROSS_BOUNDARY_SIGMA)
            M3 = M + dz_prov

        r, n_eff_grid = smoothed_residual(M3, P_g, dec, lat, lon, cont, ls_km=RESID_LS_KM)
        M_final = finalise(M, M3 + r, cont)

    # Young-age correction taper: ramp the geochem correction to zero at
    # 0 Ma (removes the present-day-bypass discontinuity for chained /
    # time-stepping consumers).  No-op for t >= YOUNG_TAPER_T_FULL, and the
    # 0 Ma slice already returned via the bypass above.  Applied before the
    # diagnostics so the reported bias/RMS/p99/Δ_rms and the written grid
    # describe the same (tapered) field.
    _r_taper = young_taper_ramp(t)
    if _r_taper < 1.0:
        M_final = apply_young_taper(t, M, M_final, cont)
        print(f"             young-age taper: ramp({t:.0f} Ma)={_r_taper:.3f} "
              f"→ correction scaled to {100*_r_taper:.0f} %")

    elapsed = time.time() - t0
    # Diagnostics
    bias_before = bias_after = rms_before = rms_after = np.nan
    if not dec.empty:
        iy = dec["iy"].astype(int).values; ix = dec["ix"].astype(int).values
        rb = dec["z"].values - M[iy, ix]; ra = dec["z"].values - M_final[iy, ix]
        bias_before = float(np.mean(rb)); bias_after = float(np.mean(ra))
        rms_before = float(np.sqrt(np.mean(rb**2))); rms_after = float(np.sqrt(np.mean(ra**2)))
    p99_before = float(np.nanpercentile(M[cont & (M>0)], 99)) if cont.any() else np.nan
    p99_after = float(np.nanpercentile(M_final[cont & (M_final>0)], 99)) if cont.any() else np.nan
    delta_rms = float(np.sqrt(np.nanmean((M_final-M)[cont]**2))) if cont.any() else np.nan

    summary = dict(
        t_Ma=int(round(t)), wall_s=elapsed,
        n_samples=int(len(df)), n_decluster=int(len(dec)) if not dec.empty else 0,
        continent_cells=int(cont.sum()),
        bias_before_m=bias_before, bias_after_m=bias_after,
        rms_before_m=rms_before,   rms_after_m=rms_after,
        p99_before_m=p99_before,   p99_after_m=p99_after,
        delta_rms_m=delta_rms,
    )
    print(f"             bias {bias_before:.0f} → {bias_after:.0f}; "
          f"RMS {rms_before:.0f} → {rms_after:.0f}; "
          f"p99 {p99_before:.0f} → {p99_after:.0f}; "
          f"Δ_rms {delta_rms:.0f}; {elapsed:.1f}s")

    if save_nc:
        out = OUT_NC_DIR / f"{int(round(t))}Ma_corrected_SW.nc"
        if out.exists():
            out.unlink()
        with nc.Dataset(out, "w") as f:
            f.createDimension("lat", len(lat)); f.createDimension("lon", len(lon))
            f.createVariable("lat", "f4", ("lat",))[:] = lat
            f.createVariable("lon", "f4", ("lon",))[:] = lon
            f.createVariable("M_orig", "f4", ("lat","lon"))[:] = M
            f.createVariable("M_corrected", "f4", ("lat","lon"))[:] = M_final
            f.createVariable("delta", "f4", ("lat","lon"))[:] = M_final - M
            f.createVariable("province", "i1", ("lat","lon"))[:] = P_g
            f.createVariable("continent_mask", "i1", ("lat","lon"))[:] = cont.astype(np.int8)
            f.createVariable("n_eff", "f4", ("lat","lon"))[:] = n_eff_grid
            f.target_age_Ma = float(t)
            f.dataset = "Scotese & Wright 2018 PaleoDEMs"
            f.config = (f"hybrid; ls={RESID_LS_KM}km; cap={DZ_MAX}m; "
                        f"n_min_p={N_MIN_P}; cob_proxy={COB_PROXY_M}m; "
                        f"local_support_km={LOCAL_SUPPORT_KM}; "
                        f"depth_fade_m={DEPTH_FADE_M}; "
                        f"classifier={PROVINCE_CLASSIFIER}")
            f.young_taper_t_full_Ma = float(YOUNG_TAPER_T_FULL)
            f.young_taper_ramp = float(young_taper_ramp(t))
    return summary, M, M_final, P_g, cont, dec


def main():
    p = argparse.ArgumentParser()
    p.add_argument("ages", nargs="*", type=float)
    p.add_argument("--all", action="store_true",
                   help="process every age available in the S&W DEM set")
    p.add_argument("--force", action="store_true",
                   help="re-run slices even if output NC already exists")
    p.add_argument("--summary-out", default=None,
                   help="path for the per-slice summary CSV")
    args = p.parse_args()

    if args.all:
        ages = sw_available_ages()
    else:
        ages = args.ages or [50.0]

    summaries = []
    skipped = 0
    for t in ages:
        out_nc = OUT_NC_DIR / f"{int(round(t))}Ma_corrected_SW.nc"
        if out_nc.exists() and not args.force:
            skipped += 1
            continue
        s, *_ = assimilate_one(float(t))
        summaries.append(s)
        if summaries:
            csv_path = (args.summary_out or
                        str(CORRECTED_DIR / "corrected_grids_SW_summary.csv"))
            pd.DataFrame(summaries).to_csv(csv_path, index=False)
    print(f"\nDone. {len(summaries)} slices processed, "
          f"{skipped} skipped (already existed).")
    print(f"NetCDFs: {CORRECTED_DIR}")


if __name__ == "__main__":
    main()
