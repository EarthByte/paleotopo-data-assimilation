# Assimilating geochemically-derived paleo-elevations into the Scotese & Wright PaleoDEM series

**Author:** prepared for Dietmar Müller
**Date:** 2026-05-10
**Inputs:**
- Scotese & Wright (2018) PaleoDEMs at 1° resolution, 0–540 Ma at 5 Myr cadence (NetCDF, variable `palaeogeography` in metres). Grid axes are paleo lat/lon — i.e. positions are in the time-slice's own reference frame, consistent with the Scotese & Wright (2018) rotation model.
- Scotese & Wright (2018) plate model — rotation file, topology / plate-boundary / plate-polygon GPML files, plus the present-day continental polygons used for sample plate-ID assignment.  All part of the same PALEOMAP PaleoAtlas release.
- Global geochemical compilation with derived crustal thickness and isostatic elevation (`Global_crustal_thickness_with_paleo_coords.csv`, ~66 600 rows in the Phanerozoic subset distributed with this release; covers 0–540 Ma). Six elevation models (Herz_0.08, Herz_0.23, Herz_0.38, Brown 2022, Davis 2009, Condie 2016) each provide low / mid / high envelope elevations driven by `z0_min/mid/max` crustal-thickness uncertainty. Reconstructed paleo-coordinates `rlat/rlon` are in the same frame as the grids.

The objective is to produce a corrected paleo-elevation field at every Myr that honours the geochemical constraints where they exist, preserves the kinematic model elsewhere, and remains geologically plausible — i.e. correction fields look like real continental relief, not noise.

---

## 1. Why this is hard

Three problems make a naive "spread the residual" interpolation unsafe:

1. **Massive temporal heterogeneity.** Sample density goes from ~33 000 in 0–50 Ma to ~200–400 in 450–540 Ma. A fixed temporal window guarantees either over-smoothing the Cenozoic or starving the Cambrian–Ordovician end of the Phanerozoic.
2. **Heavy spatial clustering.** Many samples come from the same outcrop, formation, or arc segment. Treating each row as independent over-weights well-studied terranes (e.g. SW USA, Tibet) and can drag the correction surface unphysically.
3. **Geological-shape preservation.** The kinematic model carries information about *where* mountains, cratons and basins sit; we must not destroy that pattern by smearing point corrections across province boundaries. A plateau is a flat top with steep flanks; a craton is a broad gentle bulge — the corrections must inherit these signatures.

The chosen strategy — **province-wise hypsometric rescaling with a province-aware smoothed residual field, capped and quality-weighted** — addresses all three.

---

## 2. Conceptual approach

The kinematic S&W18 map has the right *shape* — mountain belts in the right places with the right footprint — but the wrong *amplitude distribution* (its 99th-percentile elevation on land is ~2 km against an Earth-like ~4 km). The geochem data has the right *amplitude distribution* (because it's calibrated against Earth's modern crustal thickness ↔ elevation relation) but lacks spatial coverage and detail.

We therefore use the geochem data to **rescale the elevation distribution of each tectonic province on the kinematic map** so that its hypsometry matches the geochem hypsometry for the same province at the same time. The rescaling is applied as a **smooth correction field** whose spatial structure is set by where the geochem observations actually live (cluster centroids), with a province-specific correlation length.

Five elements of the design:

a. **Adaptive temporal window.** Window half-width grows with data sparsity, so the geochem hypsometry per province is always estimated from ≥ N_min samples.

b. **Cluster-aware declustering.** Samples are pre-aggregated by paleo-coord cell (and, optionally, by Tecto_Prov within the cell) so that one densely-studied outcrop becomes one observation rather than 50.

c. **Ensemble-mean elevation with explicit uncertainty.** Per-sample elevation = mean over the six geochem models; per-sample σ = max(model spread, z0_min/max envelope, Missing_Risk_DeltaT scaled). σ feeds the observation weight.

d. **Province-wise CDF matching, not point interpolation.** Within each province (Orogen, Shield, Continental Arc, Extended Crust, Platform, …) we map the kinematic-prior CDF onto the geochem CDF. This guarantees that the *distribution* of corrected elevations per province matches the data, while the *spatial pattern* of relative high/low elevations within the province is inherited from the kinematic model.

e. **Constrained smoothing.** The point-cluster residuals (geochem − kinematic) are interpolated with a province-aware Gaussian process / radial-basis kernel: short correlation length in orogens (≈300 km), long in shields (≈1000 km), and the kernel is masked to avoid bleeding across province boundaries or the COB.

---

## 3. Algorithm (per time slice t)

### Stage A — Data preparation

A1. Load `<age>_pgeog_<*>.nc` → 2-D array `M[181,361]` in metres.
A2. Load plate model rotations and COB polygons; reconstruct COB to time *t* and rasterise onto the 1° grid → continental mask `C[180,360]`.
A3. From the geochem CSV, select rows with `recon_age_Ma ∈ [t − Δt(t), t + Δt(t)]` and non-null `rlat, rlon`.

The **adaptive window** `Δt(t)` is chosen so the bin contains at least *N_min* = 200 declustered observations, capped at ±50 Myr. A reasonable default schedule (refined empirically):

| Age range | Δt (Myr) | typical sample count |
|-----------|---------:|---------------------:|
| 0–200 Ma | 5 | thousands |
| 200–500 Ma | 10 | 1k–4k |
| 500–540 Ma | 20 | 0.5k–1.5k |

### Stage B — Per-sample elevation estimate and uncertainty

For each retained row *i*:

- Mean elevation `z̄_i = mean({Herz_0.08_mid, Herz_0.23_mid, Herz_0.38_mid, Brown_mid, Davis_mid, Condie_mid})`.
- Uncertainty `σ_i² = σ_model² + σ_thickness² + σ_quality²` where
  - `σ_model = std` across the six mid models;
  - `σ_thickness = (z̄_high − z̄_low)/4` (from z0_min/mid/max envelope, treating the envelope as ~2σ);
  - `σ_quality = κ · Missing_Risk_DeltaT_km` (κ ≈ 0.3 km/km, calibrated on samples flagged `Missing_Over40pct`).
- Sample weight `w_i = 1/σ_i²`, then divided by cluster size in B1.

### B1. Spatial declustering

Bin samples into 1° (rlat, rlon) cells × Tecto_Prov. Within each (cell, province) bucket, replace the samples with a single representative:
- weighted-median elevation,
- effective σ = max(intra-cell std, mean σ_i / √n_eff),
- effective weight = sum of weights ÷ √n_eff (Kish-style penalty for clustering).

Down-weight `Missing_Over40pct = TRUE` samples by 0.25 before declustering. Samples with `Missing_Risk_DeltaT > 5 km` are dropped.

### Stage C — Province masks on the grid

We need a Tecto_Prov label for every continental grid cell. Two complementary sources:

C1. **Direct-from-data labels** where geochem samples exist: assign the dominant province of nearby samples by weighted plurality within a province-specific neighbourhood (smaller for narrower provinces — 200 km for Continental Arc, 600 km for Shield).

C2. **Geomorphic surrogate** elsewhere on the continental mask: use the prior elevation and its local 5-cell standard deviation to classify cells into {Orogen, Continental Arc, Platform/Shield, Extended Crust} via threshold rules calibrated to make the present-day map agree with C1 labels (e.g. orogen: elev > 1200 m or local σ_z > 250 m; shield: elev 0–600 m and local σ_z < 100 m).

The two are merged with C1 having priority. The province grid is a discrete categorical field `P[180,360]` defined only on the COB-mask.

### Stage D — Province-wise hypsometric rescaling

For each province *p*:

D1. Build observed CDF `F_obs,p` from the declustered geochem `(z̄, w)` pairs in *p* (within the temporal window).
D2. Build modelled CDF `F_M,p` from kinematic-prior elevations on grid cells where `P = p`.
D3. Rescaling map `T_p(z) = F_obs,p⁻¹(F_M,p(z))`.
D4. For each grid cell *c* in province *p*, compute `M'(c) = T_p(M(c))`.

This step makes the province hypsometry match the data. The relative ordering of cells inside the province is preserved — a higher prior cell remains a higher post-rescaling cell.

#### Guard rails on the rescaling map

- `T_p` is constructed only if *p* has ≥ N_min,p declustered samples (e.g. 30); otherwise `T_p = identity` (no province-wide correction; the residual interpolation in Stage E still applies locally).
- Capped: `|T_p(z) − z| ≤ Δ_max` with Δ_max = 2 km (user prior).
- Linearly tapered to identity in the lowermost / uppermost 5 % of `F_M,p` to avoid extrapolation artefacts at the tails.

### Stage E — Local residual field with province-aware smoothing

The CDF-matched field `M'` has the right hypsometry but ignores the local sample density. We refine it with a smooth correction field built from sample-level residuals.

E1. At each cluster centroid *c_j*, compute the residual `r_j = z̄_j − M'(c_j)`, with weight `w_j` and uncertainty `σ_j` from Stage B.
E2. Interpolate `r` onto the grid using a province-aware Gaussian kernel:
   `r̂(c) = Σ_j K(c, c_j) · w_j · r_j  /  Σ_j K(c, c_j) · w_j`
   with
   `K(c, c_j) = exp(−d(c,c_j)² / (2 ℓ(P(c))²))` and
   `K = 0` if `P(c) ≠ P(c_j)` (don't bleed between provinces) or if `c` is outside the COB.
   Length scales `ℓ`:
   - Continental Arc: 250 km
   - Orogen: 350 km
   - Continental Margin: 400 km
   - Extended Crust: 500 km
   - Platform: 800 km
   - Shield: 1000 km
   These are deliberately aligned with the natural along-strike / across-strike scales of those features.
E3. Effective number of samples behind each cell `N_eff(c) = (Σ K w)² / Σ (K w)²`. Where `N_eff < 1`, taper the residual smoothly to zero (no observational support).

### Stage F — Final field

`M_final(c) = M'(c) + α(c) · r̂(c)`

where `α(c)` is a confidence weight in [0,1] that uses both data support (`α → 1` as `N_eff → ∞`) and a hard cap so that `|M_final − M| ≤ Δ_max = 2 km` and `M_final ∈ [−6000, +7000]` m.

Outside the COB the field is left untouched (option to extrapolate gently in Stage G if the user wants).

### Stage G — Output

For each time slice produce:
- `M_final` — corrected paleo-elevation, NetCDF, same grid as input.
- `σ_final` — posterior standard deviation, combining the smoothed `σ_obs` with a prior σ_M (defaults to 750 m for orogens / 250 m for shields).
- `N_eff` and `P` grids — diagnostic.
- A QC plate: input vs corrected vs Δ map, hypsometric curves before/after per province, sample-residual scatter.

Stage G also generates a global `metadata.json` recording all parameters used so runs are reproducible.

---

## 4. Geological priors hard-coded

| Prior | Implementation |
|------|----------------|
| Corrections only on continent | Mask by COB at time *t*; oceanic cells of the kinematic field pass through untouched. |
| Province-specific correlation length | `ℓ(P)` table in Stage E. |
| Cap on |Δz| | `Δ_max = 2 km`, hard floor/ceiling at −6 km / +7 km. |
| Down-weight low-quality samples | `Missing_Over40pct → ×0.25`; `Missing_Risk_DeltaT` enters σ_quality and clips at 5 km. |

Optional (off by default, easy to toggle):

- Sea-level zero-crossing preservation: don't push a prior land cell below sea level via correction unless the geochem samples in that province show a substantial sub-aerial-to-flooded shift. Useful to avoid spurious epicontinental seas.
- Conservation of integrated continental hypsometry within each plate (so corrections don't add a globally implausible volume of crust).

---

## 5. Validation and diagnostics

1. **Hypsometric curves before/after** per province and globally — should converge to Earth-like for the well-sampled Phanerozoic.
2. **Leave-one-out residuals** at sample sites — corrected map should reduce sample-grid mismatch but not exactly fit (the smoothing should be apparent).
3. **Temporal continuity** — for the full pipeline, plot `M_final(c, t)` vs *t* at fixed cells across orogenic events; expect smooth growth/decay rather than 1-Myr flicker. Flicker → temporal smoothing of the residual field across adjacent maps is needed (Stage F could be extended to a 4-D Gaussian).
4. **Mass balance** — total continental volume per timestep should not jump unphysically between adjacent slices.

---

## 6. Limitations and forward path

- The geochem proxies record *crustal thickness*, not topography directly; the isostatic conversion already in the CSV assumes a particular density structure. Where lithospheric mantle removal or dynamic topography dominate (e.g. plumes, slab break-off), the proxy will be biased — visible as systematic residual sign by tectonic environment.
- For pre-700 Ma slices, sample counts are low enough that Stage D rescaling will only fire for a few provinces; the rest of the field reverts to the kinematic prior.
- The current design is independent per time slice. A natural extension is 4-D smoothing (along strike of the plate-tectonic node trajectories) so that a sample at 460 Ma also informs the 450 and 470 Ma maps with a temporal kernel.
- Province masks via the geomorphic surrogate (Stage C2) could be improved by the kinematic model's own node-level tectonic-environment tags where they exist (collisions, arcs, rifts, LIPs) — a one-time pre-processing pass over the model output.

---

## 7. Parameter table (defaults used in the proof-of-concept)

| Symbol | Meaning | Default |
|--------|--------|---------|
| Δt(t) | adaptive temporal half-window | 5/10/20/30 Myr (see table) |
| N_min | min samples per temporal bin | 200 |
| N_min,p | min declustered samples per province for D | 30 |
| Δ_max | per-cell elevation cap | 2 km |
| z_min, z_max | hard floor/ceiling | −6 km, +7 km |
| ℓ(P) | length scale per province | 250–1000 km (see Stage E) |
| α decay | residual taper as N_eff → 0 | 1 − exp(−N_eff) |
| κ | quality penalty coefficient | 0.3 km/km |
| ω_quality | weight for `Missing_Over40pct = TRUE` | 0.25 |

The proof-of-concept (`assimilate_poc.py`) implements the pipeline end-to-end on the 50 Ma slice. It writes a corrected NetCDF and a multi-panel QC PNG.
