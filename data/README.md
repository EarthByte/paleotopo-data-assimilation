# Input data

The three input datasets needed to run this workflow are distributed
via a companion Zenodo deposit so that (i) the data has its own DOI and
citation, and (ii) the git repository remains small.  All four data
components are restricted to the **Phanerozoic, 0–540 Ma**:

- 109 Scotese & Wright (2018) PaleoDEMs at 5 Myr cadence;
- the Scotese & Wright (2018) plate model + present-day continental polygons;
- the geochem CSV restricted to samples with `Age_Ma ≤ 560` (69,579 of the
  101,874 rows of the global compilation; 69,201 are 540 Ma or younger, and the
  rest are within the ±20 Myr window of the oldest slices);
- the 109 corrected NetCDFs produced by running the workflow on the
  three inputs above.

> **Zenodo DOI:** `TBA` (to be filled in once the deposit is published).

One further dataset is small enough to ship with the repository rather than
through the Zenodo deposit: the inferred craton outlines of Shirmard et al.
(2025), in `Shirmard_cratons/` (232 KB). They are cited to their own archive,
doi:10.5281/zenodo.15003746; see `Shirmard_cratons/README.md`. The independent
paleoelevation records against which the corrected grids are tested also ship
with the repository, in `independent_records/` (see its README).

After downloading the Zenodo archive, extract it into this `data/`
folder so the layout becomes:

```
data/
├── README.md                                    ← you are here
├── Scotese_Wright_2018_PaleoDEMs/
│   ├── Map01_PALEOMAP_1deg_Holocene_0Ma.nc
│   ├── Map03_PALEOMAP_1deg_Pliocene_5Ma.nc
│   ├── ...
│   └── Map88_..._540Ma.nc                       (109 NetCDFs, 5 Myr cadence; the age is read from the file name)
├── Scotese_Wright_plate_model/
│   ├── Scotese_Wright_PlateModel_2023.rot
│   ├── Scotese_Wright_PlateBoundaries_2023.gpml
│   ├── Scotese_Wright_PlatePolygons_2023.gpml
│   └── Scotese_Wright_PresentDay_ContinentalPolygons.gpml
├── geochem/
│   ├── Global_crustal_thickness_with_paleo_coords.csv
│   └── README.md                                ← geochem-dataset documentation
├── Shirmard_cratons/                            ← shipped with the repository
├── independent_records/                         ← shipped with the repository
└── corrected/                                   ← the corrected grids (also rebuilt by the workflow)
    ├── 0Ma_corrected_SW.nc
    ├── ...
    └── 540Ma_corrected_SW.nc
```

If you've placed the data elsewhere, point the workflow at it by editing
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

This takes ~10 minutes on a modern laptop for the full 109-slice sweep (the
first run also assigns plate IDs to the samples and caches them in
`geochem/sample_plate_ids_SW.npy`). A pre-built copy of `data/corrected/` is
included in the Zenodo deposit for users who want the assimilated product
directly. Each grid records its configuration in the global attribute `config`
(for the distributed grids: `resid_reach_km=250-500; delta_sigma=1.5/0;
cross_boundary_sigma=1.5; sea_sited_samples=dropped`, support radius 400 km,
craton samples excluded, declustered samples on the prior's sea excluded).
Rebuilding writes over `data/corrected/` only with `--force`; set
`PALEO_CORRECTED_DIR` to write a rebuild elsewhere and compare.
