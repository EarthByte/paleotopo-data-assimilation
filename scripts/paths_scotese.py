"""
=============================================================================
paths_scotese.py  —  Shared path configuration (Scotese & Wright workflow)
=============================================================================

Edit this one file (or set the matching environment variables) to relocate
any of the I/O destinations.  The other Scotese scripts import these
constants instead of hardcoding paths.

Inputs:
  - DEM_DIR     : Scotese & Wright 2018 PaleoDEM NetCDFs (5 Myr cadence)
  - PLATE_DIR   : Scotese & Wright revised plate model files
                  (Scotese_Wright_PlateModel_2023.rot,
                   Scotese_Wright_PlateBoundaries_2023.gpml,
                   Scotese_Wright_PlatePolygons_2023.gpml)
  - CSV_PATH    : geochem CSV with paleo-coords (rlat, rlon, Age_Ma, ...)

Outputs:
  - CORRECTED_DIR : per-Ma corrected NetCDFs + summary CSVs
  - OUTPUT_DIR    : parent of figures/ and videos/ subdirectories
                    (kept for backward compatibility; new scripts should
                    write to FIGURES_DIR or VIDEOS_DIR directly)
  - FIGURES_DIR   : all static outputs (PNG, PDF, CSV summary tables,
                    paper figures)
  - VIDEOS_DIR    : all MP4 animations
  - PAPER_FIG_DIR : the paper-numbered figures (Figures/)
  - DOCS_DIR      : Scotese-specific docs

Default project root is the parent of scripts/.  Override with
the env var PALEOTOPO_PROJECT_ROOT if you've moved the project.
=============================================================================
"""
from __future__ import annotations
import os
from pathlib import Path

# Project root: <PROJECT_ROOT>/scripts/paths_scotese.py
PROJECT_ROOT = Path(os.environ.get(
    "PALEOTOPO_PROJECT_ROOT",
    str(Path(__file__).resolve().parent.parent),
))

# Sanity fallback when running inside a session VM whose mount path differs.
if not PROJECT_ROOT.exists():
    _candidates = list(Path("/sessions").glob("*/mnt/Paleotopo_data_assimilation"))
    if _candidates:
        PROJECT_ROOT = _candidates[0]

# ---------- inputs ----------
# The layout of the Zenodo archive, extracted into data/ (see data/README.md). The second name of each pair is the
# layout of the original distributions, read when the first is absent.
def _first_existing(*cands):
    return next((c for c in cands if c.exists()), cands[0])


DEM_DIR   = _first_existing(PROJECT_ROOT / "data" / "Scotese_Wright_2018_PaleoDEMs",
                            PROJECT_ROOT / "data" / "Scotese_Wright_2018_Maps_1-88_1degX1deg_PaleoDEMS_nc_v2")
PLATE_DIR = _first_existing(PROJECT_ROOT / "data" / "Scotese_Wright_plate_model",
                            PROJECT_ROOT / "data" / "Scotese_Wright_plate_model_revised")
CSV_PATH  = _first_existing(PROJECT_ROOT / "data" / "geochem" / "Global_crustal_thickness_with_paleo_coords.csv",
                            PROJECT_ROOT / "data" / "inputs" / "crustal_thickness_topography" /
                            "Global_crustal_thickness_with_paleo_coords.csv")

# S&W-specific plate model file names
ROT_FILE        = PLATE_DIR / "Scotese_Wright_PlateModel_2023.rot"
BOUNDARIES_FILE = PLATE_DIR / "Scotese_Wright_PlateBoundaries_2023.gpml"
POLYGONS_FILE   = PLATE_DIR / "Scotese_Wright_PlatePolygons_2023.gpml"

# Continental polygons (Scotese present-day continental footprint), used
# with the S&W rotation model to produce a proper COB raster at any age.
# Searches in this order:
#   1. Scotese_Wright_PresentDay_ContinentalPolygons.gpml (current rename)
#   2. Scotese_2008_PresentDay_ContinentalPolygons.gpml
#   3. Scotese_2008_PresentDay_ContinentalPolygons.shp (legacy shapefile)
_CONT_CANDIDATES = [
    PLATE_DIR / "Scotese_Wright_PresentDay_ContinentalPolygons.gpml",
    PLATE_DIR / "Scotese_2008_PresentDay_ContinentalPolygons.gpml",
    PLATE_DIR / "Scotese_2008_PresentDay_ContinentalPolygons.shp",
]
CONTINENTAL_POLYGONS_FILE = next(
    (p for p in _CONT_CANDIDATES if p.exists()),
    _CONT_CANDIDATES[0],  # default — used to print a clear missing-file error
)

# ---------- outputs ----------
# Each can be moved with an environment variable, so a test build can be written beside the
# distributed grids rather than over them.
CORRECTED_DIR = Path(os.environ.get("PALEO_CORRECTED_DIR", str(PROJECT_ROOT / "data" / "corrected")))
OUTPUT_DIR    = Path(os.environ.get("PALEO_OUTPUT_DIR", str(PROJECT_ROOT / "outputs")))
FIGURES_DIR   = Path(os.environ.get("PALEO_FIGURES_DIR", str(OUTPUT_DIR / "figures")))
VIDEOS_DIR    = Path(os.environ.get("PALEO_VIDEOS_DIR", str(OUTPUT_DIR / "videos")))
PAPER_FIG_DIR = Path(os.environ.get("PALEO_PAPER_FIG_DIR", str(PROJECT_ROOT / "Figures")))   # the paper-numbered figures
DOCS_DIR      = PROJECT_ROOT / "docs"


def ensure_output_dirs():
    """Create writable output directories if missing."""
    for d in (CORRECTED_DIR, OUTPUT_DIR, FIGURES_DIR, VIDEOS_DIR, PAPER_FIG_DIR):
        d.mkdir(parents=True, exist_ok=True)


if __name__ == "__main__":
    print("PROJECT_ROOT     =", PROJECT_ROOT)
    for name, p in [("DEM_DIR", DEM_DIR), ("PLATE_DIR", PLATE_DIR),
                    ("ROT_FILE", ROT_FILE),
                    ("BOUNDARIES_FILE", BOUNDARIES_FILE),
                    ("POLYGONS_FILE", POLYGONS_FILE),
                    ("CONTINENTAL_POLYGONS_FILE", CONTINENTAL_POLYGONS_FILE),
                    ("CSV_PATH", CSV_PATH),
                    ("CORRECTED_DIR", CORRECTED_DIR),
                    ("OUTPUT_DIR", OUTPUT_DIR),
                    ("FIGURES_DIR", FIGURES_DIR),
                    ("VIDEOS_DIR", VIDEOS_DIR),
                    ("PAPER_FIG_DIR", PAPER_FIG_DIR),
                    ("DOCS_DIR", DOCS_DIR)]:
        print(f"{name:26s}=", p, "(exists:", p.exists(), ")")
