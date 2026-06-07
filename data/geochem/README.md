# Geochemical Dataset README

## Overview

This dataset compiles whole-rock geochemical analyses of igneous rocks alongside derived geochemical ratios, predicted crustal and mantle properties, isostatic elevation estimates, and geological/spatial metadata. It is designed to support studies of crustal thickness, mantle temperature, lithospheric structure, and tectonic provenance across different geological provinces and time periods.

Each row represents a single rock sample or averaged locality, georeferenced by latitude/longitude and assigned a tectonic province.

---

## Column Groups

### 1. Major Element Oxides (wt%)

Whole-rock major element compositions reported as weight percent oxides. Measured by XRF or EPMA.

| Column | Description |
|--------|-------------|
| `SiO2` | Silicon dioxide |
| `TiO2` | Titanium dioxide |
| `Al2O3` | Aluminium oxide |
| `FeOt` | Total iron as FeO |
| `MnO` | Manganese oxide |
| `MgO` | Magnesium oxide |
| `CaO` | Calcium oxide |
| `Na2O` | Sodium oxide |
| `K2O` | Potassium oxide |
| `P2O5` | Phosphorus pentoxide |

---

### 2. Trace Elements (ppm)

Trace element concentrations in parts per million (ppm), typically measured by ICP-MS or ICP-OES.

**Rare Earth Elements (REE):**

| Column | Element |
|--------|---------|
| `La` | Lanthanum |
| `Ce` | Cerium |
| `Pr` | Praseodymium |
| `Nd` | Neodymium |
| `Sm` | Samarium |
| `Eu` | Europium |
| `Gd` | Gadolinium |
| `Tb` | Terbium |
| `Dy` | Dysprosium |
| `Ho` | Holmium |
| `Er` | Erbium |
| `Tm` | Thulium |
| `Yb` | Ytterbium |
| `Lu` | Lutetium |

**Other Trace Elements:**

| Column | Element |
|--------|---------|
| `Sr` | Strontium |
| `Y` | Yttrium |
| `Rb` | Rubidium |
| `Ba` | Barium |
| `Hf` | Hafnium |
| `Nb` | Niobium |
| `Ta` | Tantalum |
| `Th` | Thorium |
| `U` | Uranium |
| `Pb` | Lead |
| `Zr` | Zirconium |
| `Sc` | Scandium |
| `V` | Vanadium |
| `Ni` | Nickel |
| `Cr` | Chromium |

---

### 3. Geochemical Indices and Ratios

Derived ratios used as proxies for source characteristics, melting degree, crustal thickness, and tectonic setting.

| Column | Description |
|--------|-------------|
| `A` | Aluminium saturation index or other composite index (see notes) |
| `La/Y` | Light REE enrichment vs. HREE; proxy for melt fraction and source fertility |
| `Nb/Y` | Alkalinity index; discriminates tectonic settings |
| `Ba/Sc` | Fluid-mobile to compatible element ratio; subduction proxy |
| `Nb/Yb` | Enrichment of high field strength elements relative to HREE |
| `Sr/Y` | Crustal thickness proxy; high values suggest garnet-bearing residue (thick crust or high pressure) |
| `(La/Yb)n` | Chondrite-normalised La/Yb; REE slope indicator |
| `Ce/Yb` | LREE/HREE ratio |
| `Zr/Y` | Incompatible/compatible ratio; tectonic discriminant |
| `La/Sm` | Short-range REE slope; partial melting indicator |
| `Dy/Yb` | Mid-to-heavy REE ratio; garnet vs. spinel stability indicator |
| `Sm/Yb` | Melting depth proxy |
| `Zr/Sm` | High field strength element enrichment |
| `Rb/Sr` | Crustal recycling and weathering index |
| `Nd/Yb` | REE enrichment |
| `Lu/Hf` | Radiogenic system parent/daughter ratio |
| `Ce/Y` | LREE/Y ratio |
| `Nd/Y` | MREE/Y ratio |
| `Th/Yb` | Subduction enrichment and crustal contamination proxy |
| `Gd/Yb` | Mid-REE slope; mantle source indicator |
| `Ba/V` | Fluid-mobile/compatible ratio |
| `A/CaO` | Aluminium index normalised to CaO |
| `Th/Y` | Crustal contamination indicator |
| `Ni/Sc` | Olivine vs. pyroxene fractionation proxy |
| `Cr/Sc` | Mantle source and fractionation indicator |
| `Ni/V` | Redox and fractionation proxy |
| `Cr/V` | Fractionation and source indicator |

---

### 4. Predicted Physical Properties

Model-derived estimates of crustal and lithospheric physical parameters.

| Column | Units | Description |
|--------|-------|-------------|
| `Predicted_Crustal_Thickness` | km | Crustal thickness estimated from geochemical proxies (e.g., Sr/Y, La/Yb) |
| `Mantle_Temperature_C` | °C | Estimated mantle potential temperature |
| `Crustal_Density_kg_m3` | kg/m³ | Modelled bulk crustal density |
| `Lithospheric_Mantle_Density_kg_m3` | kg/m³ | Modelled lithospheric mantle density |
| `z0_min_km` | km | Minimum estimated crustal thickness |
| `z0_mid_km` | km | Median estimated crustal thickness |
| `z0_max_km` | km | Maximum estimated crustal thickness |

---

### 5. Data Quality and Missing Value Flags

| Column | Description |
|--------|-------------|
| `Missing_Risk_DeltaT (km)` | Uncertainty in thickness estimate due to missing data (km) |
| `Missing_Rate` | Fraction of expected geochemical columns that are missing for this sample |
| `Missing_Over40pct` | Boolean flag: `TRUE` if more than 40% of geochemical columns are missing |

---

### 6. Age Information

| Column | Units | Description |
|--------|-------|-------------|
| `Age` | — | Age label or category |
| `Age_Ma` | Ma | Numeric age in millions of years before present |
| `Age error` | Ma | Analytical uncertainty on age |
| `Age_Min` | Ma | Minimum age bound |
| `Age_Max` | Ma | Maximum age bound |

---

### 7. Spatial Metadata

| Column | Description |
|--------|-------------|
| `Lat` | Latitude of sample locality (decimal degrees, WGS84) |
| `Lon` | Longitude of sample locality (decimal degrees, WGS84) |
| `Tecto_Prov` | Tectonic province classification (e.g., `Orogen`, `Shield`) |
| `Area` | Area or sub-region label (km² or descriptive, depending on context) |
| `layer` | Data layer or tectonic category |
| `continent_flag` | Integer flag indicating continental affiliation |
| `PlateID` | Tectonic plate identifier (GPlates convention) |

---

### 8. Palaeogeographic Reconstruction

Columns relating to plate tectonic reconstruction at the sample's age.

| Column | Description |
|--------|-------------|
| `recon_age_Ma` | Age used for palaeogeographic reconstruction (Ma) |
| `rlat` | Reconstructed palaeolatitude (decimal degrees) |
| `rlon` | Reconstructed palaeolongitude (decimal degrees) |

---

### 9. Isostatic Elevation Estimates

Multiple isostatic elevation models are included. Each model estimates the surface elevation (km) based on crustal and mantle density structure, using different reference mantle compositions or model parameterisations.

For each model, four columns are provided:

| Suffix | Description |
|--------|-------------|
| *(model name)* | Model identifier / crustal thickness input used |
| `*_Elev_abs_low_km` | Isostatic elevation using minimum crustal thickness (km) |
| `*_Isostatic_Elevation_absolute_km` | Isostatic elevation using median crustal thickness (km) |
| `*_Elev_abs_high_km` | Isostatic elevation using maximum crustal thickness (km) |

**Models included:**

| Model prefix | Reference |
|--------------|-----------|
| `Herz_0.23` | Herzberg et al. (mantle melt fraction = 0.23) |
| `Brown 2022` | Brown et al. (2022) |
| `Davis 2009` | Davis et al. (2009) |
| `Condie 2016` | Condie et al. (2016) |
| `Herz_0.08` | Herzberg et al. (mantle melt fraction = 0.08) |
| `Herz_0.38` | Herzberg et al. (mantle melt fraction = 0.38) |

---

## Tectonic Province Values

The `Tecto_Prov` column takes one of these nine string values
(the paleotopo-data-assimilation pipeline collapses NaN /
unrecognised values to `Other`):

- **Continental Arc** — magmatic arc built on continental crust
- **Island Arc** — magmatic arc built on oceanic crust
- **Orogen** — orogenic belt; tectonically active or recently deformed crust
- **Continental Margin** — passive or transitional continental edge
- **Extended Crust** — stretched / rifted / hyperextended continental crust
- **Basin** — sedimentary basin interior
- **Platform** — stable, gently subsiding cratonic platform
- **Shield** — stable cratonic shield; ancient, cold, thick lithosphere
- **Other** — fall-through for unclassified samples

---

## Notes

- All geochemical concentrations are in **wt%** for major oxides and **ppm** for trace elements unless otherwise stated.
- Ratios are dimensionless unless units are specified.
- `(La/Yb)n` is normalised to chondrite values (after Sun & McDonough 1989).
- Ages marked `0.00185` Ma in example rows likely represent Quaternary or very young samples.
- Rows with `Missing_Over40pct = TRUE` should be treated with caution in modelling workflows.
- Elevation values may be negative (below sea level).

---

## Phanerozoic subset for this data release

The CSV shipped in the companion Zenodo deposit is filtered to
`Age_Ma ≤ 540` — i.e. the Phanerozoic — to match the temporal scope
of the Scotese & Wright (2018) PaleoDEM series.  The full unfiltered
compilation (~99 000 continental rows, extending into the Precambrian)
was used for the Zhou et al. (2025) mohometry paper and is available
on request from the authors.

Era breakdown of the Phanerozoic subset:

| Era | Age range | Sample count |
|---|---|---|
| Cenozoic | 0–66 Ma | 32 750 |
| Mesozoic | 66–252 Ma | 22 124 |
| Paleozoic | 252–540 Ma | 11 762 |
| **Total** | **0–540 Ma** | **66 636** |

---

## How the paleotopo-data-assimilation pipeline reads this CSV

This dataset is the geochemical input to the assimilation pipeline at
https://github.com/EarthByte/paleotopo-data-assimilation.  The pipeline
consumes the following columns (the rest are kept for traceability and
downstream re-use but are not read by the assimilation step itself):

| Pipeline use | Columns |
|---|---|
| Sample age | `Age_Ma` |
| Present-day position (used for plate-ID assignment + per-slice reconstruction) | `Lat`, `Lon` |
| Tectonic-province class | `Tecto_Prov` |
| Per-sample elevation mid value (one per density model) | the six `*_Isostatic_Elevation_absolute_km` columns |
| Per-sample elevation envelope | the six `*_Elev_abs_low_km` / `_high_km` pairs |
| Crustal-thickness envelope | `z0_min_km`, `z0_mid_km`, `z0_max_km` |
| Data-quality flags | `Missing_Risk_DeltaT (km)`, `Missing_Over40pct` |

The per-sample observed elevation `z_obs` and weight `w` are computed
as (see `scripts/assimilate_scotese.py::prepare_samples`):

```
z_obs       = mean over the six "Isostatic_Elevation_absolute_km" columns
              (units: km, converted to m inside the pipeline)
σ_model     = std-dev over the same six columns
σ_thickness = inter-model RMS of (high − low) / 2
σ_quality   = 300 m/km × Missing_Risk_DeltaT (km)
σ_total     = sqrt(σ_model² + σ_thickness² + σ_quality²)
w           = 1 / σ_total²
              × 0.25 if Missing_Over40pct = True
```

The pipeline drops any row with `Missing_Risk_DeltaT > 5 km` outright.

**Note on the plate-ID / paleo-coords columns.**  The `PlateID`,
`recon_age_Ma`, `rlat`, `rlon` columns were generated against the MER21
plate-model family (Merdith et al. 2021) for the Zhou et al. (2025)
study, and are kept here for traceability.  The
paleotopo-data-assimilation pipeline does **not** use them — instead,
it assigns each sample a Scotese & Wright plate ID by point-in-polygon
test against the Scotese & Wright (2018) present-day continental polygons and
re-reconstructs the sample's `Lat`/`Lon` to each target slice age via
the Scotese & Wright (2018) rotation model.  This keeps the data
assimilation internally consistent with the Scotese & Wright kinematic
prior.

---

## Citation

If you use this dataset, please cite:

> Zhou, J., Müller, R. D. & Farahbakhsh, E. (2025).  *Machine Learning and
> Big Data Mining Reveal Earth's Deep Time Crustal Thickness and
> Tectonic Evolution: A New Chemical Mohometry Approach.*  Journal of
> Geophysical Research: Solid Earth.
> [doi:10.1029/2024JB030404](https://doi.org/10.1029/2024JB030404)

…and the Zenodo DOI of this data deposit (TBA — see the deposit page
or `../README.md`).

## License

CC-BY-4.0.  Re-use and re-distribution are permitted with attribution
to Zhou et al. (2025) and to this data deposit.  Original GEOROC and
EarthChem source data carry their own licenses; refer to those
repositories' terms if you redistribute single-source subsets.


