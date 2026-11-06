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
spatial pattern. The method builds on this paper:

Zhou, J, Farahbakhsh, E., ..., Müller, R.D., 2026. Recurrent super highlands since 2.1 Ga reveal rhythmic coupling between deep Earth and surface evolution. Geology. https://doi.org/https://doi.org/10.1130/G54718.1


---

## Highlights

- Per-slice plate-ID assignment via point-in-polygon against the Scotese and Wright
  (2018) present-day continental polygons; samples are reconstructed to
  each target slice age via the Scotese & Wright rotation model.
- The continental mask of each slice is rasterised on the sphere with pygplates'
  point-in-polygon test, so blocks that straddle the dateline or lie over a pole
  (Antarctica, Gondwana in the early Palaeozoic) are masked where they lie.
- Province-wise CDF rescaling (rank-preserving quantile mapping) +
  short-length-scale Gaussian residual kernel + continent-masked Δz
  smoothing + per-cell amplitude cap.
- Methodology includes: temporal weighting kernel, smooth
  shrinkage between N=3 and N=30 declustered samples, province pooling
  for related tectonic settings, adaptive sample-count thresholds.
- Cratons are excluded from the correction. Cratonic elevation is carried by
  sub-crustal density rather than crustal thickness (Mooney and Vidale, 2003),
  so samples inside the seismically inferred craton outlines of Shirmard et al.
  (2025) are dropped before declustering.
- Declustered observations that fall on the Scotese and Wright map's sea (prior
  elevation at or below sea level) are dropped. A sample of continental crust
  thickness says how high the land stood, not how deep the sea was; kept, these
  samples pulled the correction up over the lowlands next to them.
- The province rescaling is gated to 400 km from the nearest observation, so a
  province-wide hypsometric transfer is not applied to cells no sample can
  constrain.
- The residual correction is kept whole within 250 km of the nearest
  observation and tapered to zero at 500 km. Igneous rock chemistry samples the
  magmatic belts, not the basins between them, so the correction stays where the
  samples constrain it.
- Every corrected grid carries its support layer, `support_km`: the distance
  from each cell to the nearest declustered observation. It marks where the
  correction is whole (within 250 km), tapered (250–500 km) or absent (beyond
  500 km, where the corrected grid is the input grid).
- The grids are tested against paleoelevation estimates the correction never
  saw. In five regions (central Andes, southern Tibet, the North American
  Cordillera, the Central Pangean Mountains and eastern Asia) the input and
  corrected grids are scored against stable- and clumped-isotope paleoaltimetry,
  paleobotany, basalt vesicularity and marine datums, compiled with their DOIs
  in `data/independent_records/`. Estimates derived from igneous geochemistry
  are excluded, because igneous geochemistry is what the correction is trained on.
- 109 corrected NetCDFs for the full Phanerozoic at 5 Myr cadence,
  reproducible from a single shell command in ~5 minutes on a laptop.
- Applications: land area above 2 and 3 km through time, compared with the
  sample record of Zhou et al. (2026, Geology) on its 50 Myr bins, and the
  crustal thickness implied by the corrected elevations.

---

## Repository layout

```
paleotopo-data-assimilation/
├── README.md                      ← you are here
├── LICENSE
├── CITATION.cff                   ← citation metadata for GitHub / Zenodo
├── run_pipeline.sh                ← end-to-end runner (assimilate + diagnostics + videos)
├── build_all_figures.sh           ← rebuilds and numbers every figure and table of the paper
│
├── scripts/
│   ├── paths_scotese.py           ← shared path configuration
│   ├── sw_io.py                   ← S&W NetCDF I/O + age↔filename helpers
│   ├── plate_model_utils_scotese.py ← pyGPlates helpers (COB, SZ, provinces)
│   ├── sample_reconstruct_scotese.py ← plate-ID assignment + per-slice reconstruction
│   ├── shirmard_cratons.py        ← craton test (Shirmard et al. 2025 outlines)
│   │
│   ├── assimilate_scotese.py      ← the assimilation workflow (run this)
│   ├── build_summary_stats_scotese.py ← per-slice + global stats
│   ├── full_sweep_diagnostics_scotese.py ← Figs. 6, 7, 13
│   ├── sensitivity_refinements.py ← four-refinement ablation (Fig. S4)
│   ├── derive_crustal_thickness.py ← Airy-isostasy z_c from M_corrected (Fig. 16)
│   ├── make_sample_distribution_figure.py ← Fig. 1
│   ├── make_pipeline_illustration_figure.py ← Fig. 4 (the workflow at 100 Ma)
│   │
│   ├── esr_revision_diagnostics.py ← support layer, declustered observations, per-age areas
│   ├── esr_revision_flowchart.py  ← Fig. 2 (vector workflow diagram)
│   ├── esr_revision_maps.py       ← Fig. 5 (input, corrected, Δz with the support layer)
│   ├── esr_revision_figures.py    ← Figs. 3, 14, 15, S5, S6; Tables S2, S3
│   ├── esr_revision_regions.py    ← Figs. 8–12 and Tables S4–S8 (independent records)
│   ├── east_asia_ablation.py      ← which settings shape eastern Asia, scored on its records
│   ├── figure_checks.py           ← text-collision and label checks used by the figure scripts
│   │
│   ├── render_videos_cartopy_scotese.py ← preview MP4s (Robinson)
│   ├── render_videos_pygmt_scotese.py   ← publication MP4s (Winkel-Tripel, with the support layer)
│   └── render_video_crustal_thickness.py ← derived crustal-thickness MP4
│
├── Figures/                       ← the paper-numbered figures (Fig01..Fig16, FigS4..S6)
│   └── tables/                    ← Tables S2–S8 as CSV
│
├── docs/
│   ├── Assimilation_methodology.md ← algorithm design (math, stages, parameters)
│   └── WORKFLOW.md                 ← annotated end-to-end walkthrough
│
└── data/
    ├── README.md                  ← input data download instructions (Zenodo DOI)
    ├── Shirmard_cratons/          ← craton outlines (shipped with the repository)
    └── independent_records/       ← independent paleoelevation records (shipped with the repository)
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

### 3. Run the workflow

```bash
./run_pipeline.sh                   # ~40-50 min, full sweep + videos
./run_pipeline.sh --no-pygmt-videos # ~20 min, skips the Winkel-Tripel videos
./run_pipeline.sh --no-videos       # ~7 min, data + stats + diagnostic figures only
./run_pipeline.sh --help            # full list of options
```

The workflow is idempotent — re-running picks up where a previous run
left off.  Pass `--force` to invalidate caches and rebuild from scratch.
The distributed grids in `data/corrected/` are not overwritten without
`--force`.

### 4. Rebuild the paper's figures and tables

```bash
./build_all_figures.sh              # every figure and table (needs GMT 6.5 / pyGMT for Figs. 5, 8-12)
./build_all_figures.sh --no-pygmt   # everything except the pyGMT maps
```

The script runs each figure script in turn, stops at the first failure,
and copies the results into `Figures/` under their numbers in the paper
(see [Paper figures and tables](#paper-figures-and-tables)). Every figure
script checks its own output for overlapping text before writing it.

### 5. Outputs

| Location | Contents |
|---|---|
| `data/corrected/` | 109 corrected NetCDFs (`<age>Ma_corrected_SW.nc`) + summary CSVs and markdown dashboards |
| `outputs/` | diagnostic PNGs (`outputs/figures/`), MP4s (`outputs/videos/`), tables (`outputs/tables/`) and the support layer of every age (`outputs/esr_revision/support_km.nc`) |
| `Figures/` | the paper figures `Fig01..Fig16`, `FigS4..FigS6` and, in `Figures/tables/`, Tables S2–S8, written there by `./build_all_figures.sh` |

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
| `support_km` | f4 | great-circle distance to the nearest declustered observation (km; NaN where the age has none) |

Plus global attributes documenting the run configuration. `support_km` is
written by `assimilate_scotese.py`; for grids that do not carry it (the
v1.1.0 grids of the Zenodo archive), `scripts/esr_revision_diagnostics.py`
computes the same layer for every age into `outputs/esr_revision/support_km.nc`,
and the figure scripts read it from there.

---

## Paper figures and tables

`./build_all_figures.sh` runs the commands below in this order. To build
a single figure, run its command from `scripts/` after the first two
steps, which every later step reads.

| Paper | File in `Figures/` | Command (in `scripts/`) |
|---|---|---|
| (prerequisite) | `data/corrected/per_slice_stats_SW.csv` | `python3 build_summary_stats_scotese.py` |
| (prerequisite) | `outputs/esr_revision/{support_km.nc,declustered.csv,per_age.csv}` | `python3 esr_revision_diagnostics.py` |
| Fig. 1 | `Fig01_sample_distribution` | `python3 make_sample_distribution_figure.py` |
| Fig. 2 | `Fig02_methodology_flowchart` | `python3 esr_revision_flowchart.py` |
| Fig. 3 | `Fig03_plate_id_validation` | `python3 esr_revision_figures.py fig3` |
| Fig. 4 | `Fig04_pipeline_illustration_100Ma` | `python3 make_pipeline_illustration_figure.py` |
| Fig. 5 | `Fig05_SW_comparison` | `python3 esr_revision_maps.py fig5` |
| Figs. 6, 7, 13 | `Fig06_full_sweep_diagnostics`, `Fig07_hypsometric_curves`, `Fig13_metrics_by_era` | `python3 full_sweep_diagnostics_scotese.py` |
| Figs. 8–12, Tables S4–S8 | `Fig08_region_central_andes` … `Fig12_region_eastern_asia`, `tables/TableS4_central_andes.csv` … `tables/TableS8_eastern_asia.csv` | `python3 esr_revision_regions.py all` |
| Figs. 14, 15 | `Fig14_area_above_thresholds`, `Fig15_geology_record_comparison` | `python3 esr_revision_figures.py fig9 fig10` |
| Fig. 16 | `Fig16_derived_crustal_thickness` | `python3 derive_crustal_thickness.py --figure` |
| Fig. S4 | `FigS4_sensitivity_refinements` | `python3 sensitivity_refinements.py` |
| Figs. S5, S6, Tables S2, S3 | `FigS5_kernel_length_sweep`, `FigS6_validation_by_province`, `tables/TableS2_*`, `tables/TableS3_*` | `python3 esr_revision_figures.py figS5 tableS2` |

`esr_revision_diagnostics.py` re-runs the assimilation of every age in
memory and stops unless it reproduces the grids in `data/corrected/` to
1 m, so every derived figure belongs to the distributed grids (~2 min).
Figs. S5 and S6 and Tables S2 and S3 read the held-out cross-validation
tables (five folds of 6° blocks at every age, for kernel lengths of 100,
150, 200 and 300 km) from `outputs/kernel_length_sweep/`; the build skips
them when those tables are absent.

### The regional comparisons (Figs. 8–12, Tables S4–S8)

`esr_revision_regions.py` draws one figure per region: three ages ×
(input, corrected with the support veil and the declustered observations,
Δz on the support classes). Each independent record of
`data/independent_records/` is reconstructed to every grid age within its
age interval with the Scotese & Wright plate model, and the input and
corrected grids are read at its position (the median of the 3 × 3 cells
around it, averaged over those ages). The misfit of a grid value z
against a record that admits [lo, hi] is max(lo − z, 0, z − hi). The
record-by-record scores are written to
`outputs/tables/region_<name>_records.csv` and numbered as Tables S4–S8:

| Region | Table | Maps (Ma) | Records scored | Mean misfit, input | Mean misfit, corrected |
|---|---|---|---|---|---|
| central Andes | S4 | 50, 30, 10 | 7 | 284 m | 370 m |
| southern Tibet | S5 | 55, 35, 15 | 17 | 1,272 m | 1,122 m |
| N. American Cordillera | S6 | 95, 80, 70 | 9 | 935 m | 1,117 m |
| Central Pangean Mountains | S7 | 320, 300, 280 | 7 | 706 m | 681 m |
| eastern Asia | S8 | 160, 130, 100 | 29 | 635 m | 501 m |

These are the values of the v1.1.0 grids.
`east_asia_ablation.py` scores the eastern Asian records for each
setting of the correction in turn (`python3 east_asia_ablation.py --jobs 6`).

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

## Default workflow configuration

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
| `LOCAL_SUPPORT_KM` | 400.0 km | rescaling correction fades to zero beyond this distance from samples; override with `PALEO_LOCAL_SUPPORT_KM` |
| `PROV_TAPER_CELLS` | 2 | distance (in cells, ~110 km each) over which each province's correction fades to zero at its boundary |
| `CROSS_BOUNDARY_SIGMA` | 1.5 cells | Gaussian σ (~165 km) that bridges the post-taper neutral band at province boundaries |
| `PROVINCE_CLASSIFIER` | `"geomorphic"` | `"geomorphic"` (9 cat.) or `"hasterok"` (15 cat., Hasterok et al. 2022 ESR shapefile) |

### Observation selection

| Parameter | Default | Role |
|---|---|---|
| `DROP_CRATON_SAMPLES` | True | drop samples inside the Shirmard et al. (2025) craton outlines before declustering; override with `PALEO_DROP_CRATON` |
| `DROP_SEA_SITED_SAMPLES` | True | drop declustered samples whose cell the prior puts at or below sea level (`M_orig <= 0`); override with `PALEO_DROP_SEA_SITED` |

### Residual kernel

| Parameter | Default | Role |
|---|---|---|
| `RESID_LS_KM` | 150.0 km | per-sample Gaussian length scale for the smoothed-residual deposition |
| `RESID_REACH_KM` | (250, 500) km | the residual correction is kept whole within 250 km of the nearest declustered sample and tapered linearly to zero at 500 km; `(0, 0)` disables the taper |
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
| `DELTA_SMOOTH_SIGMA_FINAL` | 0.0 (off) | optional second Gaussian smoothing pass on Δz |
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

| Era | bias before → after (m) | p99 land before → after (m) | declustered samples, median | offshore drop, median |
|---|---|---|---|---|
| Cenozoic (5–66 Ma) | 1077 → 293 | 2360 → 3000 | 160 | 0.8 % |
| Mesozoic (66–252 Ma) | 1255 → 251 | 1840 → 2440 | 100 | 1.2 % |
| Palaeozoic (252–540 Ma) | 1332 → 513 | 2080 → 2360 | 23 | 7.5 % |

Medians over the 5 Myr slices of each era (0 Ma is the present-day bypass and
is left out); the bias is the mean of observed minus gridded elevation over the
declustered samples each slice uses. The assimilation reduces the median
per-slice bias by 750–1080 m and lifts the 99th-percentile land elevation from
a median of 2000 m to 2400 m over 5–540 Ma. Palaeozoic slices keep few
declustered samples (median 23; 13–41 in the Cambrian), and the correction there
is confined to their surroundings by the 400 km and 500 km reach limits.

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
