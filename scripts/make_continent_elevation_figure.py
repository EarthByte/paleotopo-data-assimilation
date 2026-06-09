#!/usr/bin/env python3
"""
=============================================================================
make_continent_elevation_figure.py  —  Per-continent mean elevation summary:
                                       corrected vs corrected + dyntopo diff
=============================================================================

WHAT THIS DOES
    For each age that has a combined NetCDF in --combined-dir, reconstructs
    the Scotese 2023 present-day continental polygons to that age via the
    Scotese rotation model, groups the reconstructed polygons by which
    PRESENT-DAY continent each one belongs to (using a centroid-vs-bounding-
    box check), rasterises the per-continent union onto the 1° corrected
    grid, and computes the mean elevation of M_corrected and M_combined
    within each per-continent mask.

    Output is a single xy plot with one color per continent, two line styles
    per continent:
        solid line   — M_combined   (with dyntopo time-difference correction)
        dashed line  — M_corrected  (no dyntopo correction; baseline)
    Inset legend lists the continents.  An accompanying CSV records the
    same numbers slice-by-slice for traceability.

WHY GROUP BY PRESENT-DAY CONTINENT
    The polygons are reconstructed at the target age, so they sit at their
    paleomag positions at that age and line up with M_corrected / M_combined
    cell-for-cell.  But to label each polygon with a recognisable continent
    name across deep time, we use its PRESENT-DAY centroid (via the
    polygon's get_feature().get_geometry()) and check which modern
    continental bounding box it falls into.  Fragments that have moved
    significantly through time (e.g. India before its Eocene collision)
    still get attributed to the modern continent they are now part of.

OUTPUT
    paths_scotese.OUTPUT_DIR /
        SW_continent_elevation_corrected_vs_dyntopo.png
        SW_continent_elevation_corrected_vs_dyntopo.pdf
        SW_continent_elevation_corrected_vs_dyntopo.csv

USAGE
    cd <project>/scripts_Scotese
    # Combined NetCDFs must exist first (build_dyntopo_diff_correction.py)
    python make_continent_elevation_figure.py
    python make_continent_elevation_figure.py --ages 50 100 200 300

OPTIONS
    --ages INT ...              ages to evaluate (default: all available)
    --combined-dir PATH         directory of
                                <age>Ma_corrected_plus_dyntopo_diff_young_SW.nc
    --no-pdf / --no-png / --no-csv
    --dpi INT                   PNG dpi (default 200)
    --figsize FLOAT FLOAT       figure size in inches (default 8 5)

DEPENDENCIES
    matplotlib, numpy, netCDF4, pygplates
=============================================================================
"""
from __future__ import annotations
import argparse, csv, re, sys
from pathlib import Path
import numpy as np
import netCDF4 as nc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from paths_scotese import OUTPUT_DIR, PROJECT_ROOT, CONTINENTAL_POLYGONS_FILE


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
# Modern continental bounding boxes: (lon_min, lon_max, lat_min, lat_max).
# These are deliberately broad — we're using them only to assign each
# reconstructed polygon to a modern continent label, not to constrain
# extent.  The per-age continental mask is built from the actual
# reconstructed polygon footprint at age t.
CONTINENT_BBOXES = {
    "Africa":        (-19.0,  52.0, -36.0,  38.0),
    "Antarctica":    (-180.0, 180.0, -90.0, -60.0),
    "Asia":          ( 40.0, 180.0,   5.0,  78.0),
    "Australia":     (110.0, 160.0, -44.0, -10.0),
    "Europe":        (-12.0,  40.0,  35.0,  72.0),
    "North America": (-170.0, -50.0,   8.0,  78.0),
    "South America": (-82.0, -33.0, -57.0,  13.0),
}

# Colorblind-friendlier palette; matches matplotlib's default tab10 family
# but rearranged so visually similar pairs (e.g. red/orange) don't fall on
# geographically adjacent continents.
CONTINENT_COLORS = {
    "Africa":         "#1f77b4",   # blue
    "Antarctica":     "#7f7f7f",   # grey
    "Asia":           "#ff7f0e",   # orange
    "Australia":      "#2ca02c",   # green
    "Europe":         "#9467bd",   # purple
    "North America":  "#d62728",   # red
    "South America":  "#8c564b",   # brown
}

COMBINED_FNAME_FMT = "{age}Ma_corrected_plus_dyntopo_diff_young_SW.nc"
COMBINED_FNAME_RE  = re.compile(r"^(\d+)Ma_corrected_plus_dyntopo_diff_young_SW\.nc$")
DEFAULT_COMBINED_DIR = PROJECT_ROOT / "data" / "corrected_Scotese_plus_dyntopo_diff_young"


# ---------------------------------------------------------------------------
# Polygon reconstruction + per-continent grouping
# ---------------------------------------------------------------------------
def _polygon_to_latlon(poly_on_sphere) -> np.ndarray:
    """Convert a pygplates PolygonOnSphere to an (N, 2) [lat, lon] array."""
    return np.array([p.to_lat_lon() for p in poly_on_sphere])


def _split_dateline(latlon: np.ndarray) -> list[np.ndarray]:
    """Split a polygon at the dateline if it crosses ±180° in longitude."""
    if len(latlon) < 3:
        return []
    lon = latlon[:, 1]
    dlon = np.diff(lon)
    if np.max(np.abs(dlon)) < 180:
        return [latlon]
    # Try shifting to a 0-360 longitude range; if the jump goes away, use that
    lon360 = np.where(lon < 0, lon + 360, lon)
    if np.max(np.abs(np.diff(lon360))) < 180:
        return [np.column_stack([latlon[:, 0], lon360 - 180])]
    # Otherwise split at the largest jump
    idx = int(np.argmax(np.abs(dlon))) + 1
    a, b = latlon[:idx], latlon[idx:]
    return [pts for pts in (a, b) if len(pts) >= 3]


def _present_day_centroid(feat) -> tuple[float, float] | None:
    """Return (lat, lon) of the present-day centroid of the polygon feature.

    Uses the first PolygonOnSphere returned by feat.get_geometries()."""
    for geom in feat.get_geometries():
        try:
            import pygplates as pg
        except ImportError:
            return None
        if isinstance(geom, pg.PolygonOnSphere):
            latlon = _polygon_to_latlon(geom)
            # mean-of-vertices centroid is fine for assignment purposes
            return float(latlon[:, 0].mean()), float(latlon[:, 1].mean())
    return None


def _assign_continent(centroid: tuple[float, float]) -> str | None:
    """Return the continent whose bbox contains the given (lat, lon), or None."""
    if centroid is None:
        return None
    clat, clon = centroid
    for name, (lon_min, lon_max, lat_min, lat_max) in CONTINENT_BBOXES.items():
        if lon_min <= clon <= lon_max and lat_min <= clat <= lat_max:
            return name
    return None


def reconstruct_polygons_grouped(age_ma: float):
    """Reconstruct continental polygons to age and group by present-day continent.

    Returns dict[continent_name] -> list of (N, 2) [lat, lon] arrays.
    """
    import pygplates as pg
    from plate_model_utils_scotese import get_rotation_model
    rot = get_rotation_model()
    out = []
    pg.reconstruct(str(CONTINENTAL_POLYGONS_FILE), rot, out, float(age_ma))
    grouped: dict[str, list[np.ndarray]] = {n: [] for n in CONTINENT_BBOXES}
    for rec in out:
        feat = rec.get_feature()
        centroid = _present_day_centroid(feat)
        name = _assign_continent(centroid)
        if name is None:
            continue
        rec_geom = rec.get_reconstructed_geometry()
        if not isinstance(rec_geom, pg.PolygonOnSphere):
            continue
        grouped[name].append(_polygon_to_latlon(rec_geom))
    return grouped


def build_continent_masks(grouped_polys: dict, lat: np.ndarray,
                          lon: np.ndarray) -> dict:
    """For each continent, rasterise the polygon union on the (lat, lon) grid."""
    nlat, nlon = len(lat), len(lon)
    LAT2D, LON2D = np.meshgrid(lat, lon, indexing="ij")
    pts = np.column_stack([LON2D.ravel(), LAT2D.ravel()])
    masks = {}
    for name, polys in grouped_polys.items():
        mask = np.zeros((nlat, nlon), dtype=bool)
        for latlon in polys:
            for sub in _split_dateline(latlon):
                if len(sub) < 3:
                    continue
                # MplPath wants (x, y) = (lon, lat)
                path = MplPath(sub[:, [1, 0]])
                inside = path.contains_points(pts).reshape(nlat, nlon)
                mask |= inside
        masks[name] = mask
    return masks


# ---------------------------------------------------------------------------
# Loaders + stats
# ---------------------------------------------------------------------------
def load_age(combined_dir: Path, age: int):
    f = combined_dir / COMBINED_FNAME_FMT.format(age=int(age))
    if not f.exists():
        raise FileNotFoundError(
            f"{f} not found - run `python build_dyntopo_diff_correction.py "
            f"--source young --ages {int(age)}` first")
    with nc.Dataset(f) as d:
        Mc    = np.asarray(d["M_corrected"][:], dtype=float)
        Mcomb = np.asarray(d["M_combined"][:], dtype=float)
        lat   = np.asarray(d["lat"][:], dtype=float)
        lon   = np.asarray(d["lon"][:], dtype=float)
    return Mc, Mcomb, lat, lon


def mean_elevation(M: np.ndarray, mask: np.ndarray) -> float:
    """Mean elevation in the masked region.  NaN cells excluded."""
    cells = mask & np.isfinite(M)
    n = int(cells.sum())
    if n == 0:
        return float("nan")
    return float(np.mean(M[cells]))


def auto_discover_ages(combined_dir: Path) -> list[int]:
    if not combined_dir.is_dir():
        return []
    return sorted(
        int(COMBINED_FNAME_RE.match(p.name).group(1))
        for p in combined_dir.iterdir()
        if COMBINED_FNAME_RE.match(p.name)
    )


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
def make_figure(ages: list[int], combined_dir: Path,
                figsize=(8.0, 5.0),
                write_png=True, write_pdf=True, write_csv=True, dpi=200):
    print(f"computing per-continent mean elevation for {len(ages)} ages")
    rows = []
    for age in sorted(ages, reverse=True):
        Mc, Mcomb, lat, lon = load_age(combined_dir, age)
        try:
            grouped = reconstruct_polygons_grouped(age)
        except ImportError as e:
            raise SystemExit(
                f"pygplates not installed in this Python env "
                f"({e}). Per-continent masks require pygplates.")
        masks = build_continent_masks(grouped, lat, lon)
        row = {"age": age}
        for name in CONTINENT_BBOXES:
            row[f"{name}_corrected"] = mean_elevation(Mc, masks[name])
            row[f"{name}_combined"] = mean_elevation(Mcomb, masks[name])
        rows.append(row)
        # Compact per-age trace: print just the worst-shift continent for each age
        deltas = {n: row[f"{n}_combined"] - row[f"{n}_corrected"]
                  for n in CONTINENT_BBOXES
                  if np.isfinite(row[f"{n}_combined"])
                  and np.isfinite(row[f"{n}_corrected"])}
        if deltas:
            worst = max(deltas.items(), key=lambda kv: abs(kv[1]))
            print(f"  {age:4d} Ma: largest dyntopo-driven shift  "
                  f"{worst[0]:>14s}  {worst[1]:+6.0f} m")
        else:
            print(f"  {age:4d} Ma: no continental cells found in any group")

    rows.sort(key=lambda r: r["age"])
    ages_arr = np.array([r["age"] for r in rows])

    fig, ax = plt.subplots(figsize=figsize)
    for name in CONTINENT_BBOXES:
        color = CONTINENT_COLORS[name]
        e_corr = np.array([r[f"{name}_corrected"] for r in rows])
        e_comb = np.array([r[f"{name}_combined"]  for r in rows])
        # dashed: corrected baseline
        ax.plot(ages_arr, e_corr, color=color, lw=1.2, ls="--", alpha=0.55)
        # solid: corrected + dyntopo diff
        ax.plot(ages_arr, e_comb, color=color, lw=2.0, ls="-",
                label=name)
    ax.invert_xaxis()
    ax.axhline(0, color="black", lw=0.5, ls=":", alpha=0.5)
    ax.set_xlabel("Age (Ma)", fontsize=12)
    ax.set_ylabel("Mean continental elevation (m)", fontsize=12)
    ax.grid(True, alpha=0.3, linestyle=":")
    leg = ax.legend(loc="upper right", frameon=True, fontsize=9,
                    title="Continent\n(solid = corrected + dyntopo,\n"
                          " dashed = corrected only)",
                    title_fontsize=8)
    leg.get_frame().set_alpha(0.9)
    ax.tick_params(labelsize=10)
    fig.tight_layout()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    out_base = OUTPUT_DIR / "SW_continent_elevation_corrected_vs_dyntopo"
    if write_png:
        out_png = out_base.with_suffix(".png")
        fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
        print(f"  wrote {out_png}")
    if write_pdf:
        out_pdf = out_base.with_suffix(".pdf")
        fig.savefig(out_pdf, bbox_inches="tight")
        print(f"  wrote {out_pdf}")
    plt.close(fig)
    if write_csv:
        out_csv = out_base.with_suffix(".csv")
        with open(out_csv, "w", newline="") as fh:
            w = csv.writer(fh)
            hdr = ["age_ma"]
            for name in CONTINENT_BBOXES:
                key = name.replace(" ", "_")
                hdr.append(f"{key}_mean_corrected_m")
                hdr.append(f"{key}_mean_combined_m")
            w.writerow(hdr)
            for r in rows:
                line = [r["age"]]
                for name in CONTINENT_BBOXES:
                    line.append(f"{r[f'{name}_corrected']:.1f}")
                    line.append(f"{r[f'{name}_combined']:.1f}")
                w.writerow(line)
        print(f"  wrote {out_csv}")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ages", type=int, nargs="+", default=None,
                   help="ages to evaluate. Default: ALL available in --combined-dir.")
    p.add_argument("--combined-dir", type=Path, default=DEFAULT_COMBINED_DIR,
                   help="directory of <age>Ma_corrected_plus_dyntopo_diff_young_SW.nc")
    p.add_argument("--no-png", action="store_true")
    p.add_argument("--no-pdf", action="store_true")
    p.add_argument("--no-csv", action="store_true")
    p.add_argument("--dpi", type=int, default=200)
    p.add_argument("--figsize", type=float, nargs=2, default=(8.0, 5.0),
                   metavar=("WIDTH", "HEIGHT"))
    args = p.parse_args()

    if args.ages is None:
        args.ages = auto_discover_ages(args.combined_dir)
        if not args.ages:
            raise SystemExit(
                f"ERROR: no combined NetCDFs in {args.combined_dir}. "
                f"Run `python build_dyntopo_diff_correction.py "
                f"--source young --ages <...>` first.")
        print(f"auto-discovered {len(args.ages)} ages: "
              f"{args.ages[0]}..{args.ages[-1]} Ma")
    else:
        missing = [a for a in args.ages
                   if not (args.combined_dir
                           / COMBINED_FNAME_FMT.format(age=a)).exists()]
        if missing:
            raise SystemExit(f"ERROR: no combined NetCDF for ages {missing}")

    make_figure(args.ages, args.combined_dir,
                figsize=tuple(args.figsize),
                write_png=not args.no_png,
                write_pdf=not args.no_pdf,
                write_csv=not args.no_csv,
                dpi=args.dpi)


if __name__ == "__main__":
    main()
