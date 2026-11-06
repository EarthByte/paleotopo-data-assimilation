# Independent paleoelevation records

Compilations used to test the corrected grids against paleoelevation estimates that the correction never saw. The
correction is trained on whole-rock igneous geochemistry, so only methods independent of it are admitted:
- stable-isotope (δ18O, δD) and clumped-isotope (Δ47) paleoaltimetry;
- paleobotany (CLAMP, leaf physiognomy, nearest living relative);
- basalt vesicularity;
- paleosol and fossil indicators;
- marine or coastal facies, as sea-level datums.

Crustal-thickness proxies (La/Yb, Sr/Y, Eu/Eu*) and other estimates derived from igneous geochemistry are excluded.

| file | region and interval | rows |
|---|---|---|
| `east_asia_paleoelevation_200-66Ma.csv` | eastern Asia, Jurassic–Cretaceous (and Cretaceous Tibet, `use = tibet`) | 34 |
| `central_andes_paleoelevation_60-0Ma.csv` | central Andes (Altiplano–Puna, Eastern and Western Cordillera) | 12 |
| `tibet_paleoelevation_60-0Ma.csv` | Lhasa block, Himalaya, Qiangtang, Hoh Xil | 21 |
| `cordillera_paleoelevation_100-66Ma.csv` | North American Cordillera and Western Interior Seaway margin | 9 |
| `central_pangean_mountains_paleoelevation_330-270Ma.csv` | Variscan–Appalachian belt and its forelands | 7 |
| `geology_record_comparison_50Myr_bins.csv` | the Geology record (Zhou et al., 2026) on its 50 Myr bins | 11 |

Each regional file has a `_summary.md` with what the records show and which rows are most reliable.

## Schema of the regional files

`site, region, lat, lon, coord_quality, age_min_ma, age_max_ma, elev_m, elev_min_m, elev_max_m, method,
material_or_proxy, independent_of_igneous_geochem, reference_short, doi, notes, use, lo_m, hi_m`

- **Position.** `lat`, `lon` are present-day decimal degrees; `coord_quality` is `exact` (given in the source) or `approx` (basin centroid).
- **Elevation.** `elev_m` is the published preferred value or the midpoint of the published range (the notes say which); `elev_min_m` and `elev_max_m` are the published range.
- **Interval scored.** `lo_m` and `hi_m` give the interval a record admits: the published range, `inf` or `-inf` for an open bound, and −250 to +250 m for a sea-level datum. The misfit of a grid value z is max(lo − z, 0, z − hi).
- **`use`.**
  - `strong` covers clumped isotopes, vesicularity, paleobotany with a stated error, marine datums, and δ18O with a stated lapse treatment and error.
  - `weak` covers δ18O without error, compilations, model-based or tectonic estimates, and contested readings.
  - In the East Asia file, `tibet` marks the Cretaceous Tibetan rows (reported separately). `added` marks three rows added after the scoring rule was fixed; two of them infer paleoelevation from exhumation and are reported only.
- **`doi`** was checked against CrossRef, except where the notes say otherwise.

The records are scored by `scripts/east_asia_ablation.py` (eastern Asia) and `scripts/esr_revision_regions.py` (all regions).

## The Geology record bins

`geology_record_comparison_50Myr_bins.csv` has one row per 50 Myr bin (`bin_lo`, `bin_hi`, Ma):
- `geo_p95_km`, `geo_median_km`, `geo_n`: the record's 95th percentile and median paleoelevation, and its number of samples.
- `superhighland`: whether the record identifies the bin as a super-highland interval.
- `v110_p99`, `v110_p95`, `v110_frac_gt3`: the 99th and 95th percentile land elevation and the share of land above 3 km (%) of the corrected grids, averaged over the grid ages in [bin_lo, bin_hi) without 0 Ma.
- `pre_*`: the same quantities for the correction as first submitted.
