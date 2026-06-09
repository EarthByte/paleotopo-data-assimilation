"""
=============================================================================
add_dyntopo_to_corrected_scotese.py  —  SUPERSEDED
=============================================================================

This script previously composed the geochem-corrected Scotese & Wright
paleotopography with the Young 2022 dynamic-topography anomaly
relative to present.  It is now SUPERSEDED by

    scripts_Scotese/build_dyntopo_diff_correction.py

which implements the same composition correctly: the time-difference
dyntopo(t) - dyntopo(0) is evaluated in PLATE reference frame (each
grid cell rigidly attached to its continent across time), then
cookie-cut by Scotese 2023 continental polygons and rotated into the
time-t Scotese paleomag frame via gplately.Raster.reconstruct, before
being added to M_corrected.

WHY THIS SCRIPT IS WRONG
    This script performed the dyntopo(t) - dyntopo(0) subtraction in
    the SCOTESE PALEOMAG REFERENCE FRAME (i.e. after the dyntopo had
    already been rotated by the Scotese rotation chain).  In any
    non-plate-fixed frame, continents move across time, so subtracting
    two grids at the same (lat, lon) compares dyntopo under DIFFERENT
    parts of DIFFERENT continents at the two times — physically
    meaningless.  The plate-frame subtraction in
    build_dyntopo_diff_correction.py is the correct treatment.

WHAT TO RUN INSTEAD
    cd <project>/scripts_Scotese
    python build_dyntopo_diff_correction.py --source young \\
        --ages 5 25 50 75 100 150 200 250 300

Output directory and NetCDF schema differ from this old script:
    OLD: data/corrected_Scotese_plus_dyntopo/<age>Ma_corrected_plus_dyntopo_SW.nc
         variables: M_corrected, z_dyntopo, z_dyntopo_anomaly, M_combined, continent_mask
    NEW: data/corrected_Scotese_plus_dyntopo_diff_young/<age>Ma_corrected_plus_dyntopo_diff_young_SW.nc
         variables: M_corrected, z_dyntopo_diff, M_combined  (CF-1.10 + ACDD-1.3)

The figure scripts (make_comparison_figures_dyntopo.py,
make_flooded_fraction_figure.py, make_dyntopo_panels_figure.py) have
all been retargeted at the new directory and variable names; running
this stub will simply error out with a clear pointer to the
replacement.
=============================================================================
"""
from __future__ import annotations
import sys


def main() -> int:
    print(__doc__, file=sys.stderr)
    print("\nERROR: This script is superseded by "
          "scripts_Scotese/build_dyntopo_diff_correction.py.\n"
          "       Run that script instead — see the header docstring "
          "above for the exact command.",
          file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main())
