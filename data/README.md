# Input data

The three input datasets needed to run this pipeline are distributed
via a companion Zenodo deposit so that (i) the data has its own DOI and
citation, and (ii) the git repository remains small.  All four data
components are restricted to the **Phanerozoic, 0–540 Ma**:

- 109 Scotese & Wright (2018) PaleoDEMs at 5 Myr cadence;
- the Scotese & Wright (2018) plate model + present-day continental polygons;
- the geochem CSV filtered to samples with `Age_Ma ≤ 540` (~66 600
  samples out of the full ~99 000-row global compilation);
- the 109 corrected NetCDFs produced by running the pipeline on the
  three inputs above.

> **Zenodo DOI:** `TBA` (to be filled in once the deposit is published).

After downloading the Zenodo archive, extract it into this `data/`
folder so the layout becomes:

```
data/
├── README.md                                    ← you are here
├── Scotese_Wright_2018_PaleoDEMs/
│   ├── 0_paleodem.nc
│   ├── 5_paleodem.nc
│   ├── ...
│   └── 540_paleodem.nc                          (109 NetCDFs, 5 Myr cadence)
├── Scotese_Wright_plate_model/
│   ├── Scotese_Wright_PlateModel_2023.rot
│   ├── Scotese_Wright_PlateBoundaries_2023.gpml
│   ├── Scotese_Wright_PlatePolygons_2023.gpml
│   └── Scotese_Wright_PresentDay_ContinentalPolygons.gpml
├── geochem/
│   ├── Global_crustal_thickness_with_paleo_coords.csv
│   └── README.md                                ← geochem-dataset documentation
└── corrected/                                   ← populated by the pipeline
    ├── 0Ma_corrected_SW.nc
    ├── ...
    └── 540Ma_corrected_SW.nc
```

If you've placed the data elsewhere, point the pipeline at it by editing
`scripts/paths_scotese.py` or by setting the environment variable
`PALEOTOPO_PROJECT_ROOT` to the parent folder.

## Citation and licensing of the inputs

The three input products carry their own citations and licenses:

| Product | Citation | License |
|---|---|---|
| Scotese & Wright (2018) PALEOMAP PaleoAtlas — PaleoDEMs + rotation model + plate-boundary / plate-polygon / continental-polygon GPML files | Scotese, C. R. & Wright, N. M. (2018). PALEOMAP PaleoAtlas for GPlates and the PaleoData Plotter Program. PALEOMAP Project. | Creative Commons (consult PALEOMAP for current terms) |
| Whole-rock geochemical compilation | Zhou, J., Müller, R. D. & Farahbakhsh, E. (2025). *JGR: Solid Earth*. [doi:10.1029/2024JB030404](https://doi.org/10.1029/2024JB030404) | See `geochem/README.md` for terms |

You must cite each of these in any work derived from this repository.
The Apache 2.0 license that governs this codebase does **not** override
the licenses of the input data.

## Reproducing the corrected NetCDFs

Once the inputs are in place, the corrected NetCDFs are produced by:

```bash
cd scripts
python assimilate_scotese.py --all
```

This takes ~5 minutes on a modern laptop for the full 109-slice sweep.
A pre-built copy of `data/corrected/` is also included in the Zenodo
deposit for users who want the assimilated product directly without
re-running the pipeline.
