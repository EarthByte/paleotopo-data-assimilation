"""
=============================================================================
make_supermountains_figure.py  —  Figure 8 of the Earth-Science Reviews paper
=============================================================================

Computes the Phanerozoic "supermountains index" from the corrected S&W
NetCDFs and plots a time-series figure suitable for the paper.

DEFINITION
    S(t) = (area of continental cells with M_corrected > THRESHOLD_M) /
           (area of continental cells at time t)

    Area is computed via cos(lat) weighting on the 1° × 1° grid.  The
    threshold defaults to 3000 m (high mountains: Andes / Tibet / Cordillera);
    a secondary series at 2000 m is also plotted as a faint reference.

OUTPUT
    paper/Scotese/Fig08_supermountains_index.png

ANNOTATIONS
    Three classical supermountain epochs labelled on the time series:
      - Late Cenozoic (Himalayan-Tibetan + Andean + Cordilleran)
      - Late Permian-Triassic (Variscan-Appalachian / Hercynian, central Pangaea)
      - Cambrian-Ordovician (Pan-African / Trans-Gondwanan)

USAGE
    cd <project>/scripts_Scotese
    python make_supermountains_figure.py
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

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from paths_scotese import CORRECTED_DIR
PROJECT_ROOT = HERE.parent
FIG_DIR = PROJECT_ROOT / "paper" / "Scotese"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Configuration
THRESHOLD_PRIMARY = 3000.0   # "supermountains"
THRESHOLD_SECONDARY = 2000.0 # high-elevation reference band

# Geological era boundaries (Ma) for shading
ERAS = [
    (0,   66,  "Cenozoic",   "#fae6c8"),
    (66,  252, "Mesozoic",   "#cfe6cf"),
    (252, 540, "Paleozoic",  "#cfd9e6"),
]

# Supermountain epoch labels (centre Ma, label, label-y-offset relative to series max)
EPOCHS = [
    ( 25, "Cordilleran /\nAlpine-Himalayan\n(Cenozoic)"),
    (270, "Variscan-Appalachian\n/ Hercynian"),
    (505, "Pan-African /\nTrans-Gondwanan"),
]


def cell_area_weights(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    """Per-cell area weights (sphere) for a regular lat-lon grid.
    Normalised so that ocean+continent total area = 1.0."""
    LAT2D, _ = np.meshgrid(lat, lon, indexing="ij")
    w = np.cos(np.radians(LAT2D))
    return w / w.sum()


def supermountain_fraction(M_corr, cont, lat, lon, threshold_m: float) -> float:
    w = cell_area_weights(lat, lon)
    cont_area = (cont * w).sum()
    if cont_area <= 0:
        return np.nan
    high = (cont & (M_corr > threshold_m)) * w
    return float(high.sum() / cont_area)


def build_timeseries(threshold_m: float):
    ages, frac, cont_area = [], [], []
    files = sorted(CORRECTED_DIR.glob("*Ma_corrected_SW.nc"),
                   key=lambda f: int(f.name.split("Ma")[0]))
    for f in files:
        with nc.Dataset(f) as d:
            try: t = int(round(float(d.target_age_Ma)))
            except AttributeError: t = int(f.name.split("Ma")[0])
            Mc = d.variables["M_corrected"][:].astype(float)
            cont = d.variables["continent_mask"][:].astype(bool)
            lat = d.variables["lat"][:].astype(float)
            lon = d.variables["lon"][:].astype(float)
        ages.append(t)
        frac.append(supermountain_fraction(Mc, cont, lat, lon, threshold_m))
        w = cell_area_weights(lat, lon)
        cont_area.append(float((cont * w).sum()))
    df = pd.DataFrame({"t_Ma": ages, "S_t": frac, "cont_area_frac": cont_area})
    return df.sort_values("t_Ma").reset_index(drop=True)


def main():
    print("Building supermountains time series …")
    df_3km = build_timeseries(THRESHOLD_PRIMARY)
    df_2km = build_timeseries(THRESHOLD_SECONDARY)
    out_csv = FIG_DIR.parent.parent / "data" / "corrected_Scotese" / "supermountains_index.csv"
    df_3km.merge(df_2km.rename(columns={"S_t": "S_t_2km"})[["t_Ma","S_t_2km"]],
                 on="t_Ma").to_csv(out_csv, index=False)
    print(f"  wrote {out_csv}")

    fig, ax = plt.subplots(figsize=(11, 5))

    # shaded era backgrounds
    ymax = max(np.nanmax(df_3km.S_t * 100), np.nanmax(df_2km.S_t_2km * 100 if "S_t_2km" in df_2km else df_2km.S_t * 100))
    ymax = max(ymax * 1.20, 4)
    for lo, hi, name, col in ERAS:
        ax.axvspan(lo, hi, color=col, alpha=0.55, zorder=0)
        ax.text((lo+hi)/2, 1.02, name, transform=ax.get_xaxis_transform(),
                ha="center", va="bottom", fontsize=9, color="#333")

    # supermountain time series (2 km = thin grey reference, 3 km = thick red)
    ax.fill_between(df_2km.t_Ma, df_2km.S_t * 100, 0,
                    color="#888888", alpha=0.25,
                    label="continental area > 2000 m  (reference)")
    ax.plot(df_3km.t_Ma, df_3km.S_t * 100,
            color="#c0392b", lw=2.2,
            label="continental area > 3000 m  (supermountains)")

    # epoch labels — place at series-local maxima
    for epoch_ma, label in EPOCHS:
        sub = df_3km[(df_3km.t_Ma >= epoch_ma - 25) & (df_3km.t_Ma <= epoch_ma + 25)]
        if sub.empty:
            continue
        local_max = sub.S_t.max() * 100
        ax.annotate(label,
                    xy=(epoch_ma, local_max),
                    xytext=(epoch_ma, max(local_max + 1.5, 3.0)),
                    ha="center", va="bottom",
                    fontsize=9, color="#333",
                    arrowprops=dict(arrowstyle="-", color="#666", lw=0.6))

    ax.set_xlabel("Age (Ma)", fontsize=11)
    ax.set_ylabel("Continental area fraction (%)", fontsize=11)
    ax.set_title("Phanerozoic supermountains index from the corrected paleotopography grids",
                 fontsize=12, pad=14)
    ax.invert_xaxis()
    ax.set_xlim(540, 0)
    ax.set_ylim(0, ymax)
    ax.grid(True, alpha=0.3, lw=0.6)
    ax.legend(loc="upper left", framealpha=0.9, fontsize=9)
    plt.tight_layout()

    out_png = FIG_DIR / "Fig08_supermountains_index.png"
    plt.savefig(out_png, dpi=200, bbox_inches="tight")
    print(f"  wrote {out_png}")


if __name__ == "__main__":
    main()
