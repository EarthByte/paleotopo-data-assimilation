"""
=============================================================================
derive_crustal_thickness.py
=============================================================================

Inverts Airy isostasy on the assimilated paleo-elevation field to estimate
crustal thickness z_c (km) at every continental cell.  Uses standard
single-pair densities (no per-cell density model) — the result is
suitable for first-order continental-scale visualisations and as an
input field to other isostatically-driven models, NOT as a calibrated
mohometry replacement.

For a continental column in Airy compensation with reference column at
sea level (h = 0, z_c = z_0) and compensation depth below the deepest
crustal root:

    h ≥ 0  (subaerial):
        z_c = z_0 + h × ρ_m / (ρ_m − ρ_c)

    h < 0  (flooded continental cell, water of density ρ_w):
        z_c = z_0 + h × (ρ_m − ρ_w) / (ρ_m − ρ_c)

In both cases z_c grows with h, but the subaerial formula uses the
full mantle density contrast whereas the sub-aqueous one is damped by
the water on top of the depression.

Default constants (overrideable at the top of the file or via CLI):

    z_0  = 33.0  km      normal continental crust at sea level
    ρ_c  = 2850. kg/m³   bulk continental crust
    ρ_m  = 3300. kg/m³   asthenospheric mantle
    ρ_w  = 1030. kg/m³   sea water

These are the textbook values used in most paleo-elevation /
paleo-crustal-thickness conversions.  They are deliberately *spatially
uniform*; the geochem-mohometry inputs to the assimilation use
per-sample density models (Herzberg, Brown, Davis, Condie), but here we
operate on the gridded corrected elevation field where no cell-level
density information is available.

USAGE
    cd <repo root>/scripts
    python derive_crustal_thickness.py 100               # one slice
    python derive_crustal_thickness.py --all             # all 109 slices
    python derive_crustal_thickness.py --figure          # also build the
                                                         # 6-panel preview
                                                         # figure (50, 100,
                                                         # 200, 300, 400,
                                                         # 500 Ma)

OUTPUTS
    data/corrected/<age>Ma_crustal_thickness_SW.nc
        lat, lon  (f4)   paleo-coordinate axes
        z_c       (f4)   crustal thickness in km (NaN over ocean)
        h_m       (f4)   M_corrected echoed back (m) for traceability
        continent_mask   (i1)
    paper/Scotese/Fig11_derived_crustal_thickness.png   (--figure only)
=============================================================================
"""
from __future__ import annotations
import argparse, sys
from pathlib import Path
import numpy as np
import netCDF4 as nc

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from paths_scotese import CORRECTED_DIR, OUTPUT_DIR

# ---------------------------------------------------------------------------
# Standard Airy-isostasy constants — edit here if you want to swap
# density structures.
# ---------------------------------------------------------------------------
Z0_KM    = 33.0     # normal continental crust at sea level
RHO_C    = 2850.0   # bulk continental crust  (kg/m³)
RHO_M    = 3300.0   # asthenospheric mantle   (kg/m³)
RHO_W    = 1030.0   # sea water               (kg/m³)

# Physical clipping of the output.  Crustal thickness < 0 is unphysical,
# > ~80 km is rare (Himalaya / Tibet upper bound) — both extremes here
# come from continental cells that sit at the M_corrected hard floor
# (−11 km) / ceiling (+9 km).
Z_C_MIN_KM = 0.0
Z_C_MAX_KM = 80.0


def airy_thickness(h_m, continent_mask):
    """Airy crustal-thickness inversion on a 2-D elevation field.

    Parameters
    ----------
    h_m : array, m
        Surface elevation in metres (signed; negative = below sea level).
    continent_mask : array, bool
        Polygon-derived continental footprint at this age.  We **do
        not** use this as the visualisation mask, because the Scotese
        2008 polygons reconstructed via S&W23 rotations diverge from
        the S&W18 PaleoDEM's own paleo-coastline at deep time (most
        visibly before ~350 Ma), creating "ghost continent" shapes in
        the rendered thickness map.  Instead the output is masked to
        the actual paleo-land of the rendered elevation grid
        (h_m ≥ 0).  Submerged continental shelves are excluded — z_c
        is reported only where the cell is subaerial in the elevation
        field that drives the inversion.  `continent_mask` is accepted
        for API back-compat but only used to disambiguate at the COB
        edges (a cell is land if h ≥ 0 AND the polygon agrees, or
        h ≥ 0 alone).

    Returns
    -------
    z_c_km : array, km
        Crustal thickness in km, with NaN over ocean cells and over
        any cells that aren't paleo-subaerial in `h_m`.
    """
    h_km = h_m / 1000.0
    subaerial = (h_km >= 0.0)
    z_c = np.full_like(h_km, np.nan, dtype=np.float64)

    # h >= 0: full mantle contrast
    z_c[subaerial] = Z0_KM + h_km[subaerial] * RHO_M / (RHO_M - RHO_C)

    # Clip to a physically plausible range and mask off everything
    # below sea level (including submerged continental shelves — this
    # is the change that suppresses the polygon-derived ghost-continent
    # shadows in the rendered thickness videos).
    z_c = np.clip(z_c, Z_C_MIN_KM, Z_C_MAX_KM)
    z_c = np.where(subaerial, z_c, np.nan)
    return z_c.astype(np.float32)


# ---------------------------------------------------------------------------
# Per-slice driver
# ---------------------------------------------------------------------------
def process_one(age_ma: int):
    src = CORRECTED_DIR / f"{int(age_ma)}Ma_corrected_SW.nc"
    if not src.exists():
        raise FileNotFoundError(src)
    with nc.Dataset(src) as d:
        lat  = np.asarray(d.variables["lat"][:],            dtype=np.float32)
        lon  = np.asarray(d.variables["lon"][:],            dtype=np.float32)
        Mc   = np.asarray(d.variables["M_corrected"][:],    dtype=np.float32)
        cont = np.asarray(d.variables["continent_mask"][:], dtype=bool)

    z_c = airy_thickness(Mc, cont)

    out = CORRECTED_DIR / f"{int(age_ma)}Ma_crustal_thickness_SW.nc"
    with nc.Dataset(out, "w") as d:
        d.createDimension("lat", lat.size)
        d.createDimension("lon", lon.size)
        d.createVariable("lat", "f4", ("lat",))[:] = lat
        d.createVariable("lon", "f4", ("lon",))[:] = lon
        v = d.createVariable("z_c", "f4", ("lat", "lon"))
        v.units = "km"
        v.long_name = "crustal thickness, Airy isostasy on M_corrected"
        v[:] = z_c
        v2 = d.createVariable("h_m", "f4", ("lat", "lon"))
        v2.units = "m"
        v2.long_name = "corrected paleo-elevation (echo of M_corrected)"
        v2[:] = Mc
        v3 = d.createVariable("continent_mask", "i1", ("lat", "lon"))
        v3[:] = cont.astype(np.int8)
        d.target_age_Ma   = float(age_ma)
        d.z_0_km          = Z0_KM
        d.rho_c_kgm3      = RHO_C
        d.rho_m_kgm3      = RHO_M
        d.rho_w_kgm3      = RHO_W
        d.formula         = ("h>=0: z_c = z_0 + h*rho_m/(rho_m-rho_c); "
                             "h<0:  z_c = z_0 + h*(rho_m-rho_w)/(rho_m-rho_c)")
        d.history = "produced by scripts/derive_crustal_thickness.py"
    # Diagnostic
    land = z_c[np.isfinite(z_c)]
    if land.size:
        print(f"  {age_ma:>3} Ma  →  z_c  min/median/p99/max = "
              f"{land.min():.1f} / {np.median(land):.1f} / "
              f"{np.percentile(land, 99):.1f} / {land.max():.1f} km  "
              f"(n={land.size})")
    return out


# ---------------------------------------------------------------------------
# Preview figure  (loaded lazily so the data-only driver doesn't import
# matplotlib + cartopy)
# ---------------------------------------------------------------------------
def make_preview_figure(ages=(50, 100, 200, 300, 400, 500)):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    import cartopy.crs as ccrs

    proj = ccrs.Robinson()
    norm = mcolors.Normalize(vmin=20, vmax=70)
    cmap = plt.get_cmap("YlOrBr")

    fig, axes = plt.subplots(3, 2, figsize=(12, 10),
                             subplot_kw={"projection": proj})
    pcm = None
    for ax, t in zip(axes.flat, ages):
        f = CORRECTED_DIR / f"{int(t)}Ma_crustal_thickness_SW.nc"
        if not f.exists():
            process_one(t)
        with nc.Dataset(f) as d:
            lat = np.asarray(d.variables["lat"][:])
            lon = np.asarray(d.variables["lon"][:])
            z_c = np.asarray(d.variables["z_c"][:], dtype=float)
            h_m = np.asarray(d.variables["h_m"][:], dtype=float)
        ax.set_global()
        LON2D, LAT2D = np.meshgrid(lon, lat)
        pcm = ax.pcolormesh(LON2D, LAT2D, z_c, cmap=cmap, norm=norm,
                            transform=ccrs.PlateCarree(),
                            shading="auto", rasterized=True)
        # Sea-level coastline contour on the elevation grid the z_c
        # field was derived from (matches the rendered topography).
        ax.contour(LON2D, LAT2D, (h_m >= 0).astype(float),
                   levels=[0.5], colors="black", linewidths=0.4,
                   transform=ccrs.PlateCarree())
        ax.set_title(f"{t} Ma", fontsize=11)
        ax.gridlines(lw=0.25, color="grey", alpha=0.4)

    cax = fig.add_axes([0.22, 0.05, 0.56, 0.018])
    cb = fig.colorbar(pcm, cax=cax, orientation="horizontal", extend="both")
    cb.set_label("Derived crustal thickness z_c  (km)", fontsize=10)
    fig.suptitle(
        f"Derived crustal thickness from Airy isostasy on corrected S&W18 elevation"
        f"\n(z₀={Z0_KM:g} km, ρ_c={RHO_C:g}, ρ_m={RHO_M:g}, ρ_w={RHO_W:g} kg/m³)",
        fontsize=12, y=0.98
    )
    fig.subplots_adjust(left=0.02, right=0.98, top=0.94, bottom=0.10,
                        hspace=0.18, wspace=0.05)

    out_dir = OUTPUT_DIR.parent / "paper" / "Scotese"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_png = out_dir / "Fig11_derived_crustal_thickness.png"
    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_png}")
    return out_png


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("age", nargs="?", type=int,
                   help="single slice age in Ma (omit + use --all for all slices)")
    p.add_argument("--all", action="store_true",
                   help="process every corrected NetCDF in data/corrected/")
    p.add_argument("--figure", action="store_true",
                   help="also build the 6-panel preview figure "
                        "(paper/Scotese/Fig11_derived_crustal_thickness.png)")
    p.add_argument("--ages", type=int, nargs="+",
                   help="figure-only: choose your own six ages "
                        "(default: 50 100 200 300 400 500)")
    args = p.parse_args()

    if args.all:
        files = sorted(CORRECTED_DIR.glob("*Ma_corrected_SW.nc"),
                       key=lambda f: int(f.name.split("Ma")[0]))
        ages = [int(f.name.split("Ma")[0]) for f in files]
        print(f"processing {len(ages)} slices …")
        for t in ages:
            process_one(t)
    elif args.age is not None:
        print(f"processing {args.age} Ma …")
        process_one(args.age)
    elif not args.figure:
        p.error("specify an age, --all, or --figure")

    if args.figure:
        make_preview_figure(tuple(args.ages or (50, 100, 200, 300, 400, 500)))


if __name__ == "__main__":
    main()
