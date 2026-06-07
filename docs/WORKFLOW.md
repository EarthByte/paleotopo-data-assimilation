# Workflow — Scotese & Wright assimilation, end-to-end

End-to-end walkthrough for the standalone S&W workflow.  Assumes:
- the repo has been cloned to `~/Documents/Software/Paleotopo_data_assimilation/`,
- the S&W PaleoDEMs and revised plate-model files are in `data/`,
- the geochem CSV is in `data/inputs/`,
- Python and the pyGPlates / cartopy / pygmt / ffmpeg dependencies are installed.

If anything below references a path or value you want to change, edit
**`paths_scotese.py`** — it's the single point of truth for I/O locations.

## 0. Verify your environment

```bash
cd ~/Documents/Software/Paleotopo_data_assimilation/scripts
python paths_scotese.py
```

Every path printed should end with `(exists: True)`.  If `CONTINENTAL_POLYGONS_FILE`
shows `False`, point the search list in `paths_scotese.py` at whatever
filename you have.

```bash
python -c "import pygplates, cartopy; print(pygplates.__version__, cartopy.__version__)"
which gmt ffmpeg
```

## 1. Conceptual overview

The pipeline operates one age-slice at a time.  For each age *t*:

```
                                            ┌─────────────────────┐
   S&W PaleoDEM NetCDF (181×361)  ────────▶│                     │
                                            │                     │
   geochem CSV (rlat/rlon paleo-coords)  ─▶│ assimilate_scotese  │
        + Tecto_Prov tags                   │                     │
                                            │ Stage A — load      │
   S&W rotations + plate topologies      ─▶│ Stage B — declust   │
   + Scotese 2008 continental polygons     │ Stage C — provs     │ ─▶ <t>Ma_corrected_SW.nc
                                            │ Stage D — CDF rescale     (M_orig, M_corrected,
                                            │ Stage E — residual         delta, province,
                                            │ Stage F — cap & finalise   continent_mask, n_eff)
                                            │                     │
                                            └─────────────────────┘
```

The scientific design (math behind each stage) is identical to the
is described in `../docs/Assimilation_methodology.md`.

## 2. Run the production assimilation

```bash
cd scripts
python assimilate_scotese.py --all
```

Per-slice runtime is 1–2 s; full S&W sweep (107 ages) ≈ 4 min.  Output
NetCDFs go to `data/corrected/<age>Ma_corrected_SW.nc`.

The pipeline is **idempotent** — re-running `--all` skips slices whose
NetCDF already exists.  Add `--force` to overwrite.

To run a subset:
```bash
python assimilate_scotese.py 50 100 200 500   # explicit list
python assimilate_scotese.py --all --force    # full sweep, overwrite
```

The script writes an incremental summary CSV next to the NetCDFs at
`data/corrected/corrected_grids_SW_summary.csv`.

## 3. Build time-dependent summary stats

```bash
python build_summary_stats_scotese.py
```

Four files appear in `data/corrected/`:

| file | what it gives you |
|------|-------------------|
| `per_slice_stats_SW.csv` | the primary time-dependent table — one row per Ma, all metrics |
| `province_delta_summary_SW.csv` | one row per (age, province) — Δz by tectonic province through time |
| `before_after_summary_SW.md` | human-readable global dashboard |
| `era_summary_SW.md` | era-binned aggregate stats (Cenozoic / Mesozoic / Paleozoic) |

This step is fast (< 1 minute).

## 4. Build temporal-evolution figures

```bash
python full_sweep_diagnostics_scotese.py
```

Three figures land in `outputs/`:

| figure | what it shows |
|--------|---------------|
| `SW_full_sweep_diagnostics.png` | 4-panel time-series: p99 land elevation, Δ RMS / max, input + corrected hypsometric heatmaps |
| `SW_hypsometry_selected_ages.png` | overlaid corrected hypsometric curves at 8 selected ages vs Earth-modern |
| `SW_metrics_by_era.png` | per-era box plots of bias, RMS, p99, Δz RMS — useful for paper figures |

## 5. Render videos

### Preview (no GMT required — Robinson projection)

```bash
python render_videos_cartopy_scotese.py --fps 12
```

Outputs in `outputs/`:
- `SW_paleotopo_elevation_540-0Ma_preview.mp4`     — corrected elevation
- `SW_paleotopo_original_540-0Ma_preview.mp4`      — original S&W input
- `SW_paleotopo_delta_540-0Ma_preview.mp4`         — Δz
- `SW_paleotopo_elevation_samples_540-0Ma_preview.mp4` — corrected + controls

### Publication-grade Winkel-Tripel (pyGMT)

```bash
python render_videos_pygmt_scotese.py --fps 10
```

Same four mode files, no `_preview` suffix, ~3 s per frame.

If you only want a subset:
```bash
python render_videos_pygmt_scotese.py --modes original elevation_samples
python render_videos_pygmt_scotese.py --ages 0 50 100 200 500
```

## 6. Time-dependent stats — recommended interpretation

For each metric in `per_slice_stats_SW.csv` here's the recommended
narrative angle:

| Metric | What it tells you |
|--------|-------------------|
| `n_decluster` | Sample-density indicator.  When this drops below ~30 the pipeline falls back to small residual-only corrections.  Phanerozoic always has plenty; data sparsity rises sharply in the Cambrian–Ordovician. |
| `bias_before_m` / `bias_after_m` | The kinematic prior's amplitude bias against the geochem evidence.  Positive ⇒ model is too low.  Reduction from before → after is the pipeline's main numerical claim. |
| `rms_before_m` / `rms_after_m` | Sample-to-grid scatter.  Always reduces less than the bias because the residual kernel deliberately doesn't try to fit each cluster to zero error — fitting noise is bad practice. |
| `land_corr_p99_m` | Top of the corrected map.  Compare against ~4 km Earth-modern reference.  Drops in data-sparse intervals are physically reasonable, not artefacts. |
| `delta_rms_m`, `delta_max_m` | Correction magnitude per slice.  Bounded by the 2 km cap.  Spikes correspond to large prior errors (often hothouse / orogen-rich intervals). |
| `era` column | Use for grouping; era-binned aggregates are usually cleaner story than time-series for paper tables. |

## 7. Customising / extending

Single-point edits to common knobs:

| Change | File |
|--------|------|
| Move the project | `paths_scotese.py` (or `PALEOTOPO_PROJECT_ROOT` env var) |
| Different cap / length scale / N_min,p | `assimilate_scotese.py` top constants |
| Different temporal window schedule | `dt_half()` in `assimilate_scotese.py` |
| Different province classification rule | `province_grid()` in `plate_model_utils_scotese.py` |
| Different COB rasterisation source | `cob_mask()` in `plate_model_utils_scotese.py` (set `use_polygons=False` to fall back to DEM threshold) |
| Drop subduction-zone overlays in videos | `DRAW_SZ = False` in `render_videos_pygmt_scotese.py` |
| Reconstruct samples via S&W rotations | see "Sample paleo-coords" below |
