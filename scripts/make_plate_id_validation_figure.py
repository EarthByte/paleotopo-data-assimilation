"""
=============================================================================
make_plate_id_validation_figure.py  —  Figure 4 of the Earth-Science Reviews paper
=============================================================================

Three-panel validation of the plate-ID assignment and per-slice
reconstruction step:

  (a) Global map of geochemical samples at present-day, coloured by their
      assigned S&W plate ID.
  (b) Sample-trajectory polylines for representative samples: a few
      "marker" samples chosen from each major continent, with their
      reconstructed Lat/Lon plotted as polylines from 540 Ma → present.
  (c) Offshore-drop time series: fraction of in-window samples whose
      reconstructed position falls outside the slice's continent mask,
      with a 25% warning threshold.

OUTPUT
    Figures/Fig04_plate_id_validation.png
=============================================================================
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import netCDF4 as nc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from paths_scotese import CSV_PATH, CORRECTED_DIR
import assimilate_scotese as A
from sample_reconstruct_scotese import ScoteseSampleReconstructor

PROJ_ROOT = HERE.parent
FIG_DIR = PROJ_ROOT / "Figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)


def main():
    print("Loading geochem CSV + assigned plate IDs …")
    df = A.get_geochem()
    print(f"  {len(df)} samples after dropping oceanic")

    fig = plt.figure(figsize=(13, 12))

    # ----------------------- (a) sample plate IDs at present-day -----------------------
    ax_a = fig.add_subplot(2, 1, 1, projection=ccrs.Robinson())
    ax_a.set_global()
    ax_a.coastlines(lw=0.4, color="grey")
    ax_a.gridlines(lw=0.3, color="grey", alpha=0.3)
    # Limit to plate IDs with enough samples to be visually distinguishable,
    # then collapse the rest to a "minor plates" bucket so the legend stays
    # short.  Each major plate gets its own colour.
    pid_counts = df["sw_plate_id"].value_counts()
    major = pid_counts.head(8).index.tolist()
    df["pid_cat"] = df["sw_plate_id"].apply(lambda p: p if p in major else -1)
    cmap = plt.cm.tab10
    colors = {pid: cmap(i) for i, pid in enumerate(major)}
    colors[-1] = (0.7, 0.7, 0.7, 0.8)
    plate_names = {
        101: "North America", 201: "South America", 301: "Eurasia",
        302: "Eastern Eurasia", 304: "Iberia", 501: "India",
        604: "South China",   611: "Indochina",   801: "Australia",
        701: "Africa",        901: "Antarctica",  -1:  "minor plates",
    }
    for pid in major + [-1]:
        sub = df[df["pid_cat"] == pid]
        if sub.empty: continue
        ax_a.scatter(sub["Lon"], sub["Lat"], c=[colors[pid]],
                     s=5, alpha=0.55, edgecolor="none",
                     transform=ccrs.PlateCarree(),
                     label=f"{int(pid) if pid != -1 else 'minor'} — {plate_names.get(pid, '?')}"
                            if pid != -1 else "minor plates")
    ax_a.set_title("(a)  Geochemical samples coloured by assigned S&W plate ID (present-day positions)",
                   fontsize=11, loc="left")
    ax_a.legend(loc="lower left", fontsize=8, ncol=2, framealpha=0.95,
                 markerscale=3)

    # ----------------------- (b) sample trajectories -----------------------
    ax_b = fig.add_subplot(2, 2, 3, projection=ccrs.Robinson())
    ax_b.set_global()
    ax_b.coastlines(lw=0.4, color="grey")
    ax_b.gridlines(lw=0.3, color="grey", alpha=0.3)
    # Pick one marker sample from each major continent so the trajectory plot is legible
    sr = ScoteseSampleReconstructor()
    markers = {
        "N. America": (40.0, -100.0, 101),
        "S. America": (-15.0, -60.0, 201),
        "Africa":     ( 10.0,  20.0, 701),
        "India":      ( 22.0,  78.0, 501),
        "Australia":  (-25.0, 135.0, 801),
        "Antarctica": (-80.0,   0.0, 901),
    }
    ages_traj = list(range(0, 541, 20))
    color_cycle = plt.cm.tab10(np.linspace(0, 1, len(markers)))
    for i, (name, (la0, lo0, pid)) in enumerate(markers.items()):
        df_pt = pd.DataFrame({"Lat":[la0], "Lon":[lo0], "sw_plate_id":[pid]})
        lats, lons = [], []
        for t in ages_traj:
            rla, rlo = sr.reconstruct(df_pt, t)
            lats.append(rla[0]); lons.append(rlo[0])
        ax_b.plot(lons, lats, "-", color=color_cycle[i], lw=1.5,
                  transform=ccrs.Geodetic(), label=name)
        ax_b.scatter([lons[0]], [lats[0]], color=color_cycle[i],
                     s=40, edgecolor="black", linewidths=0.5,
                     transform=ccrs.PlateCarree(), zorder=5)
        ax_b.text(lons[0], lats[0]+3, "0", color=color_cycle[i], fontsize=8,
                  transform=ccrs.PlateCarree())
        ax_b.text(lons[-1], lats[-1]+3, "540",
                  color=color_cycle[i], fontsize=8,
                  transform=ccrs.PlateCarree())
    ax_b.set_title("(b)  Marker-sample reconstruction trajectories 0–540 Ma\n(circle = present-day)",
                   fontsize=11, loc="left")
    ax_b.legend(loc="lower left", fontsize=8, framealpha=0.95)

    # ----------------------- (c) offshore drop fraction vs time -----------------------
    ax_c = fig.add_subplot(2, 2, 4)
    summary_csv = CORRECTED_DIR / "per_slice_stats_SW.csv"
    if summary_csv.exists():
        stats = pd.read_csv(summary_csv)
    else:
        stats = None

    # Compute offshore-drop fraction per slice on-the-fly so we don't rely
    # on the summary file (which may not record it).
    from sw_io import nearest_cell_index
    rows = []
    files = sorted(CORRECTED_DIR.glob("*Ma_corrected_SW.nc"),
                   key=lambda f: int(f.name.split("Ma")[0]))
    for f in files:
        with nc.Dataset(f) as d:
            try: t = int(round(float(d.target_age_Ma)))
            except AttributeError: t = int(f.name.split("Ma")[0])
            cont = d.variables["continent_mask"][:].astype(bool)
            lat  = d.variables["lat"][:]
            lon  = d.variables["lon"][:]
        sub_raw = A.prepare_samples(df, t)
        n_raw = len(sub_raw)
        if n_raw == 0:
            rows.append({"t_Ma": t, "n_raw": 0, "frac_off": np.nan})
            continue
        iy = nearest_cell_index(lat, sub_raw["rlat"].values)
        ix = nearest_cell_index(lon, sub_raw["rlon"].values)
        # Matches the assimilator's accounting: polygon mask only.
        on = cont[iy, ix]
        rows.append({"t_Ma": t, "n_raw": n_raw,
                     "frac_off": 100 * (n_raw - on.sum()) / n_raw})

    dfr = pd.DataFrame(rows).sort_values("t_Ma")
    ax_c.plot(dfr["t_Ma"], dfr["frac_off"], color="#c0392b", lw=1.2)
    ax_c.fill_between(dfr["t_Ma"], dfr["frac_off"], 0, color="#c0392b", alpha=0.18)
    ax_c.axhline(25, color="grey", lw=0.8, ls="--",
                 label="25 % warning threshold")
    ax_c.set_xlabel("Age (Ma)")
    ax_c.set_ylabel("offshore-drop fraction (%)")
    ax_c.set_title("(c)  Fraction of samples landing outside the continent mask after reconstruction",
                   fontsize=11, loc="left")
    ax_c.invert_xaxis()
    ax_c.grid(True, alpha=0.3, lw=0.6)
    ax_c.legend(loc="upper left", fontsize=9, framealpha=0.95)

    plt.tight_layout()
    out = FIG_DIR / "Fig04_plate_id_validation.png"
    plt.savefig(out, dpi=180, bbox_inches="tight")
    print(f"  wrote {out}")


if __name__ == "__main__":
    main()
