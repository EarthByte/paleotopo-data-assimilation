# paleotopo-data-assimilation

> Open-source data-assimilation workflow that corrects the elevation
> amplitudes of kinematic paleotopography models using a global
> compilation of geochemically-derived paleo-elevation estimates.
> Applied here to the Scotese & Wright (2018) PaleoDEM series at 5 Myr
> intervals across the Phanerozoic.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository implements the methodology described in:

> Zhou, J., Müller, R. D. & Farahbakhsh, E. (in review).  *Assimilating
> crustal-thickness-derived elevations into PaleoDEMs: an open framework
> for amplitude-corrected Phanerozoic paleotopography.*  Earth-Science
> Reviews.

The kinematic Scotese & Wright (2018) PaleoDEMs (S&W18) have a reasonable
*shape* (location of mountain belts, continental footprint) but 
underestimate the *amplitude* of past topography — their
top 1 % of land is at ~2 km versus Earth's modern ~4 km.  This workflow
ingests ~67 000 geochemically-derived paleo-elevation samples and
updates the S&W18 grids on a slice-by-slice basis, recovering an
Earth-modern-like hypsometric envelope while preserving the kinematic
spatial pattern.

---

## Highlights

- Per-slice plate-ID assignment via point-in-polygon against the Scotese and Wright
  (2018) present-day continental polygons; samples are reconstructed to
  each target slice age via the Scotese & Wright rotation model.
- Province-wise CDF rescaling (rank-preserving quantile mapping) +
  short-length-scale Gaussian residual kernel + continent-masked Δz
  smoothing + per-cell amplitude cap.
- Methodology includes: temporal weighting kernel, smooth
  shrinkage between N=3 and N=30 declustered samples, province pooling
  for related tectonic settings, adaptive sample-count thresholds.
- 109 corrected NetCDFs for the full Phanerozoic at 5 Myr cadence,
  reproducible from a single shell command in ~5 minutes on a laptop.
- Demonstration application: a spatially-resolved Phanerozoic
  *supermountains index* time series.

---

## Repository layout

```
paleotopo-data-assimilation/
├── README.md                      ← you are here
├── LICENSE                        ← MIT
├── CITATION.cff                   ← citation metadata for GitHub / Zenodo
├── run_pipeline.sh                ← end-to-end runner (assimilate + diagnostics + videos)
├── build_all_figures.sh           ← single-command rebuild of every Fig01..Fig11
│
├── scripts/
│   ├── paths_scotese.py           ← shared path configuration
│   ├── sw_io.py                   ← S&W NetCDF I/O + age↔filename helpers
│   ├── plate_model_utils_scotese.py ← pyGPlates helpers (COB, SZ, provinces)
│   ├── sample_reconstruct_scotese.py ← plate-ID assignment + per-slice reconstruction
│   │
│   ├── assimilate_scotese.py      ← PRODUCTION pipeline (run this)
│   ├── build_summary_stats_scotese.py ← per-slice + global stats
│   ├── full_sweep_diagnostics_scotese.py ← temporal-evolution figures
│   ├── sensitivity_refinements.py ← four-refinement ablation (Fig 10)
│   ├── derive_crustal_thickness.py ← Airy-isostasy z_c from M_corrected (Fig 11)
│   │
│   ├── make_*.py                  ← figure scripts (Fig 2/3/4/5/8)
│   ├── draw_methodology_flowchart.py ← Fig 1
│   │
│   ├── render_videos_cartopy_scotese.py ← preview MP4s (Robinson)
│   ├── render_videos_pygmt_scotese.py   ← publication MP4s (Winkel-Tripel)
│   └── render_video_crustal_thickness.py ← derived crustal-thickness MP4
│
├── Figures/                       ← publication-ready PNG + PDF set
│   ├── Fig01_methodology_flowchart.png
│   ├── Fig02_sample_distribution.png
│   ├── Fig03_pipeline_illustration_100Ma.png
│   ├── Fig04_plate_id_validation.png
│   ├── Fig05[abc]_SW_comparison_*.{png,pdf}
│   ├── Fig06_full_sweep_diagnostics.png
│   ├── Fig07_hypsometric_curves.png
│   ├── Fig08_supermountains_index.png
│   ├── Fig09_metrics_by_era.png
│   └── Fig10_sensitivity_refinements.png
│
├── docs/
│   ├── Assimilation_methodology.md ← algorithm design (math, stages, parameters)
│   └── WORKFLOW.md                 ← annotated end-to-end walkthrough
│
└── data/
    └── README.md                  ← input data download instructions (Zenodo DOI)
```

The `data/` subfolders for inputs and outputs are excluded from git via
`.gitignore`; see `data/README.md` for how to populate them.

---

## Quick start

### 1. Clone and set up the Python environment

```bash
git clone https://github.com/EarthByte/paleotopo-data-assimilation.git
cd paleotopo-data-assimilation

# Recommended: a clean conda env
conda create -n paleotopo python=3.11 numpy pandas scipy netCDF4 xarray \
    matplotlib cartopy pygplates pygmt -c conda-forge
conda activate paleotopo
brew install ffmpeg                # macOS  (or `conda install -c conda-forge ffmpeg`)
```

### 2. Download the input data

See [`data/README.md`](data/README.md).  Inputs (Scotese & Wright
PaleoDEMs, the S&W revised plate model, and the geochem CSV) are
distributed via Zenodo to keep the repository small and to give the data
its own citation.

### 3. Run the pipeline

```bash
./run_pipeline.sh                   # ~40-50 min, full sweep + videos
./run_pipeline.sh --no-pygmt-videos # ~20 min, skips the Winkel-Tripel videos
./run_pipeline.sh --no-videos       # ~7 min, data + stats + paper figures only
./run_pipeline.sh --help            # full list of options
```

The pipeline is idempotent — re-running picks up where a previous run
left off.  Pass `--force` to invalidate caches and rebuild from scratch.

### 4. Outputs

| Location | Contents |
|---|---|
| `data/corrected/` | 109 corrected NetCDFs (`<age>Ma_corrected_SW.nc`) + summary CSVs and markdown dashboards |
| `outputs/` | preview MP4s, publication MP4s, diagnostic PNGs |
| `Figures/` | the eleven figures `Fig01..Fig11` written there by `./build_all_figures.sh` |

Each `<age>Ma_corrected_SW.nc` contains:

| variable | dtype | description |
|---|---|---|
| `lat`, `lon` | f4 | paleo-coordinate axes (−90..+90 / −180..+180) |
| `M_orig` | f4 | original S&W18 elevation (m) |
| `M_corrected` | f4 | assimilated elevation (m) |
| `delta` | f4 | M_corrected − M_orig (m) |
| `province` | i1 | tectonic-province index (0–8; −1 = ocean) |
| `continent_mask` | i1 | continental footprint (1 = land) |
| `n_eff` | f4 | effective sample support per cell |

Plus global attributes documenting the run configuration.

---

## Default pipeline configuration

All knobs live at the top of `scripts/assimilate_scotese.py` and can be
tuned without code changes elsewhere.  The defaults below are the
production values used for the Phanerozoic sweep released with the paper:

### Per-cell elevation envelope

| Parameter | Default | Role |
|---|---|---|
| `DZ_MAX` | 2000 m | per-cell \|Δz\| cap |
| `ZMIN`, `ZMAX` | −11000, +11000 m | hard floor / ceiling for M_corrected |
| `SUBAERIAL_FLOOR_M` | 1.0 m | cells subaerial in M_orig are kept ≥ this in M_corrected |
| `DEPTH_FADE_M` | 500.0 m | corrections fade linearly to zero at this depth; smooth coastline transition |

### Per-province CDF rescaling

| Parameter | Default | Role |
|---|---|---|
| `N_MIN_P` | 3 | min Kish-effective samples per province for any rescaling |
| `N_FULL_P` | 30 | full-amplitude rescaling above this n_eff (smooth shrinkage in between) |
| `PROVINCE_POOLING` | True | allow related provinces to pool samples (down-weight 0.3) |
| `LOCAL_SUPPORT_KM` | 1000.0 km | rescaling correction fades to zero beyond this distance from samples |
| `PROV_TAPER_CELLS` | 2 | distance (in cells, ~110 km each) over which each province's correction fades to zero at its boundary |
| `CROSS_BOUNDARY_SIGMA` | 3.0 cells | Gaussian σ that bridges the post-taper neutral band at province boundaries |
| `PROVINCE_CLASSIFIER` | `"geomorphic"` | `"geomorphic"` (9 cat.) or `"hasterok"` (15 cat., Hasterok et al. 2022 ESR shapefile) |

### Residual kernel

| Parameter | Default | Role |
|---|---|---|
| `RESID_LS_KM` | 150.0 km | per-sample Gaussian length scale for the smoothed-residual deposition |
| `KAPPA_QUALITY` | 300 m/km | σ contribution from `Missing_Risk_DeltaT` |
| `W_MISSING_OVER40` | 0.25 | down-weight factor when `Missing_Over40pct = True` |

### Temporal kernel

| Parameter | Default | Role |
|---|---|---|
| `TEMPORAL_KERNEL` | `"triangular"` | sample temporal weight; `"uniform"`, `"triangular"`, or `"gaussian"` |
| `dt_half(t)` | 5 / 10 / 20 / 30 Myr | adaptive temporal half-window (0–200 / 200–500 / 500–800 / 800–1000 Ma) |

### Post-processing of the Δz field

| Parameter | Default | Role |
|---|---|---|
| `DELTA_MEDIAN_WINDOW` | 5 | 5×5 median-filter despeckle on Δz (kills clusters up to ~12 pixels) |
| `DELTA_SMOOTH_SIGMA` | 1.5 cells | first Gaussian smoothing pass on Δz (~165 km) |
| `DELTA_SMOOTH_SIGMA_FINAL` | 4.5 cells | final Gaussian smoothing pass on Δz (~495 km) |
| `FINAL_LAND_SMOOTH_SIGMA` | 0.0 (disabled) | optional subaerial-only Gaussian on M_corrected |
| `CONT_MASK_CLOSING_ITER` | 6 | morphological closing iterations on the continent mask (fills polygon gaps up to ~12 cells wide) |

### Continent / province classification

| Parameter | Default | Role |
|---|---|---|
| `COB_PROXY_M` | −1500 m | fallback DEM threshold for the continent mask if polygons are unavailable |
| `ARC_DIST_KM` | 450 km | "Continental Arc" province threshold around subduction zones (geomorphic classifier only) |
| `MARGIN_DIST_KM` | 200 km | "Continental Margin" width from the COB (geomorphic classifier only) |

To override any of these for a one-off run, edit the constants at the
top of `scripts/assimilate_scotese.py` (or `scripts/plate_model_utils_scotese.py`
for the COB / arc / margin thresholds) and re-run with `--force`.

---

## Headline results

| Era | bias before → after (m) | p99 land before → after (m) | n_decluster median | offshore drop median |
|---|---|---|---|---|
| Cenozoic (0–66 Ma) | 1722 → 640 | 2349 → 3980 | 231 | 4 % |
| Mesozoic (66–252 Ma) | 3257 → 290 | 1840 → 3618 | 198 | 8 % |
| Paleozoic (252–540 Ma) | 3732 → 800 | 2120 → 3668 | 118 | 22 % |

The assimilation reduces the median per-slice sample-residual bias by
500–1100 m and lifts the 99th-percentile land elevation from a typical
2 km to ~3.7 km across the Phanerozoic.

---

## Citation

If you use this code or the corrected grids, please cite both:

1. **The methodology paper** (this repository's reference):

   > Zhou, J., Müller, R. D. & Farahbakhsh, E. (in review). *Assimilating
   > crustal-thickness-derived elevations into PaleoDEMs: an open
   > framework for amplitude-corrected Phanerozoic paleotopography.*
   > Earth-Science Reviews.

2. **The underlying mohometry**:

   > Zhou, J., Müller, R. D. & Farahbakhsh, E. (2025). *Machine Learning
   > and Big Data Mining Reveal Earth's Deep Time Crustal Thickness and
   > Tectonic Evolution: A New Chemical Mohometry Approach.* Journal of
   > Geophysical Research: Solid Earth.
   > [doi:10.1029/2024JB030404](https://doi.org/10.1029/2024JB030404)

A `CITATION.cff` file is included for automated citation tools (GitHub's
"Cite this repository" button and Zenodo).

The underlying plate model is based on:

Scotese, C. R. & Wright, N. M. (2018).  PALEOMAP PaleoAtlas for
  GPlates and the PaleoData Plotter Program.  PALEOMAP Project.


---

## License

MIT License — see [`LICENSE`](LICENSE).

Copyright © 2026 The EarthByte Group, University of Sydney.

The third-party inputs distributed via Zenodo carry their own licenses;
see `data/README.md` for terms.
