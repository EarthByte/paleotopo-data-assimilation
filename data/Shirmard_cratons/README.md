# Shirmard et al. (2025) inferred craton outlines

`Craton_Boundaries_Inferred.shp` — 20 named craton polygons, 5.00e7 km2 total,
in EPSG:4326.

    Shirmard, H., et al. (2025). Craton boundary detection from full-waveform
    tomography model reveals links to critical metal deposits.
    Geoscience Frontiers 16, 102176. doi:10.1016/j.gsf.2025.102176

    Data archive: doi:10.5281/zenodo.15003746

These files are redistributed here under the terms of the archive's licence so
the assimilation runs without a separate download; the Zenodo archive remains
the citable source and should be cited in preference to this copy.

Used by `scripts/shirmard_cratons.py` to exclude geochemical samples falling
inside a craton before declustering. Cratonic surface elevation is carried by
sub-crustal density rather than crustal thickness (Mooney and Vidale, 2003,
Lithos 71, 185-193), so a crustal-thickness-to-elevation conversion has no skill
there. See `docs/Assimilation_methodology.md` section 2(f).

The companion `Craton_Boundaries.shp` from the same archive holds 546 detected
boundary line segments rather than closed polygons and is not used here.
