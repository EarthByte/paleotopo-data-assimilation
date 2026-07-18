# paleotopo-data-assimilation

> Open-source data-assimilation workflow that corrects the elevation
> amplitudes of kinematic paleotopography models using a global
> compilation of geochemically-derived paleo-elevation estimates.
> Applied here to the Scotese & Wright (2018) PaleoDEM series at 5 Myr
> intervals across the Phanerozoic.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository implements the methodology described in:

> Zhou, J., Müller, R. D. & Farahbakhsh, E. (in prep).  *Assimilating
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
spatial pattern. The method builds on this paper:

Zhou, J, Farahbakhsh, E., ..., Müller, R.D., 2026. Recurrent super highlands since 2.1 Ga reveal rhythmic coupling between deep Earth and surface evolution. Geology. https://doi.org/https://doi.org/10.1130/G54718.1


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

## Optional: dynamic-topography composition test

A supplementary analysis (paper Supp Mat §S.1) composes the
geochem-corrected paleotopography with the Young et al. (2022) GLD428
dynamic-topography model across 50–300 Ma.  Scripts:

| Script | Output |
|---|---|
| `scripts/build_dyntopo_diff_correction.py` | `<age>Ma_corrected_plus_dyntopo_diff_young_SW.nc` — corrected + plate-frame dyntopo time-difference rotated to Scotese paleomag frame |
| `scripts/make_dyntopo_panels_figure.py` | Supp Fig S1 — six-panel mantle-flow correction overview |
| `scripts/make_comparison_figures_dyntopo.py` | Supp Fig S2 — three Fig 5-style 2 × 3 comparisons at (300 \| 250), (200 \| 150), (100 \| 50) Ma |
| `scripts/make_flooded_fraction_figure.py` | Supp Fig S3 — continental flooding through time |
| `scripts/make_continent_elevation_figure.py` | Supp Fig S4 — per-continent mean elevation through time |
| `scripts/render_videos_pygmt_dyntopo_scotese.py` | publication-quality Winkel-Tripel MP4s; modes `combined` / `dyntopo_diff` / `dyntopo_absolute` / `corrected` |
| `scripts/render_corrected_plus_dyntopo_video.py` | cartopy / Robinson preview MP4s of the composed field |

### Data — Young et al. (2022) dynamic topography grids

Download the GLD428 grids from EarthByte (University of Sydney):

| Release | URL | Files |
|---|---|---|
| 5-Myr cadence, 0–250 Ma (Merdith2021 plate model) | [`/Dynamic_Topography/gld428/`](https://www.earthbyte.org/webdav/ftp/Dynamic_Topography/gld428/) | `gld428-PlateFrame-0-250Ma.zip` (48 MB), `gld428-MantleFrame-0-250Ma.zip` (67 MB) |
| 20-Myr cadence (Merdith2021 plate model), 0–1000 Ma | [`/Dynamic_Topography/gld428_m21/`](https://www.earthbyte.org/webdav/ftp/Dynamic_Topography/gld428_m21/) | `gld428_plate_frame_grids.zip` (88 MB), `gld428_mantle_frame_grids.zip` (46 MB) |

The dynamic-topography models can also be interactively visualised via the
GPlates Portal Dynamic Topography page (under "Young et al., 2022"):
<https://portal.gplates.org/portal/dt/>.

`build_dyntopo_diff_correction.py` reads the plate-frame grids; the
optional upstream `rotate_young_dyntopo_to_scotese.py` reads the
mantle-frame grids.  Unzip directly into the project's `data/` directory:

```
data/
├── Young2022_gld428_grids_5Myr/         # from gld428-PlateFrame-0-250Ma.zip
├── Young2022_gld428_grids_20Myr/        # from gld428_plate_frame_grids.zip
└── …
```

That layout is the default the scripts expect — no further configuration
needed.  If you'd rather keep the grids somewhere else, override at run
time:

```bash
# Option A — environment variable (points at a parent data dir):
export PALEOTOPO_DATA_ROOT=/path/to/your-data-root

# Option B — per-invocation CLI flag:
python scripts/build_dyntopo_diff_correction.py \
    --young-dir /path/to/Young2022_gld428_grids_20Myr \
    --source young --ages $(seq 0 5 300) --step-myr 5
```

Cite the data as:

> Young, A., Flament, N., Williams, S. E., Merdith, A., Cao, X., &
> Müller, R. D. (2022). Long-term Phanerozoic sea level change from
> solid Earth processes. *Earth and Planetary Science Letters*, **584**,
> 117451.

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

### Young-age correction taper

| Parameter | Default | Role |
|---|---|---|
| `YOUNG_TAPER_T_FULL` | 15 Ma | age by which the geochem correction reaches full strength; it is ramped linearly to zero at 0 Ma (`0.0` disables the taper) |

At 0 Ma the assimilation is bypassed by design — the present-day DEM is
observed directly and must not be altered, so `delta ≡ 0`.  Without a
taper the correction jumps to its full amplitude (`delta ≈ +100 m` mean
over land) at the very next slice (5 Ma), so a chained or time-stepping
consumer of the corrected DEMs (e.g. a stepwise landscape-evolution
model) sees a broad continental step appear across the 0↔5 Ma boundary
that is an artefact of the bypass rather than tectonics.  The method was
designed to amplitude-correct the kinematic prior in *deep* time, where
sparse geochem samples are spatially bled into a geologically plausible
field; the youngest Cenozoic — plates near their modern positions, S&W
already close to truth — never required it.  `YOUNG_TAPER_T_FULL` ramps
the correction `ramp(t) = min(1, t / T_full)` so it grows smoothly from
the fixed 0 Ma DEM instead of switching on abruptly.  Because it is a
per-age, same-frame rescale of `M_corrected − M_orig`, it is independent
of plate motion, leaves ages ≥ `T_full` unchanged, and only affects the
5 and 10 Ma slices at the default value.

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
