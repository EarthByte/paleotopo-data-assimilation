#!/usr/bin/env bash
# =============================================================================
# run_pipeline.sh  —  end-to-end assimilation runner
# =============================================================================
#
# Sequentially runs:
#   1. assimilate all 109 slices (0–540 Ma, 5 Myr cadence)
#   2. build summary stats CSVs + markdown dashboards
#   3. build temporal-evolution diagnostic figures (Figs 6, 7, 9)
#   4. run the four-refinement ablation (Fig 10)
#   5. derive Airy-isostasy crustal thickness from M_corrected (Fig 11)
#   6. render preview videos (cartopy, Robinson, 4 modes)
#   7. render publication videos (pygmt, Winkel-Tripel, 4 modes)
#
# Designed to be launched under nohup so it survives ssh disconnect.
# All stages are idempotent — a re-run skips slices/videos that already
# exist unless `--force` is given.
#
# USAGE
#   cd <repo root>
#   nohup ./run_pipeline.sh > pipeline.log 2>&1 &
#   disown
#   tail -f pipeline.log
#
# OPTIONS
#   --no-videos          skip ALL video rendering (cartopy AND pygmt)
#   --no-pygmt-videos    skip just the slow Winkel-Tripel pygmt videos
#   --no-sensitivity     skip the four-refinement ablation (Fig 10)
#   --no-thickness       skip the Airy-isostasy crustal-thickness derivation
#                        (Fig 11)
#   --force              re-do anything previously cached:
#                          - plate-ID caches (re-assigned from polygons)
#                          - per-slice corrected NetCDFs
#                          - video frame PNGs
#   --pygmt-cadence N    for pyGMT video runs use cadence N Ma (default 5).
#   -h, --help           show this help and exit
#
# RUNTIME (single-threaded, no `--force`, starting from clean run)
#   assimilate (109 slices)   ≈  5 min
#   stats + diagnostics       ≈  1 min
#   sensitivity ablation      ≈  1 min
#   cartopy previews          ≈ 10 min
#   pyGMT publication videos  ≈ 20–30 min (cadence 5)
#   ----------------------------------------------------------------
#   Total                     ≈ 40–50 min
#
# =============================================================================
set -u
set -o pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

NO_VIDEOS=0
NO_PYGMT=0
NO_SENSITIVITY=0
NO_THICKNESS=0
FORCE=0
PYGMT_CADENCE=5
while [[ $# -gt 0 ]]; do
    case "$1" in
        --no-videos)            NO_VIDEOS=1 ;;
        --no-pygmt-videos)      NO_PYGMT=1 ;;
        --no-sensitivity)       NO_SENSITIVITY=1 ;;
        --no-thickness)         NO_THICKNESS=1 ;;
        --force)                FORCE=1 ;;
        --pygmt-cadence)        PYGMT_CADENCE="$2"; shift ;;
        --pygmt-cadence=*)      PYGMT_CADENCE="${1#*=}" ;;
        -h|--help)
            awk '/^[^#]/{exit} NR>1' "$0" | sed 's/^# \{0,1\}//'
            exit 0
            ;;
        *)  echo "ERROR: unknown option $1 (try --help)"; exit 2 ;;
    esac
    shift
done

ts()  { date '+%Y-%m-%d %H:%M:%S'; }
log() { printf '[%s] %s\n' "$(ts)" "$*"; }
hdr() { printf '\n[%s] === %s ===\n' "$(ts)" "$*"; }
run_tagged() {
    local tag="$1"; shift
    "$@" 2>&1 | sed -u "s|^|[$tag] |"
    return "${PIPESTATUS[0]}"
}

PY=$(command -v python3 || command -v python)
if [[ -z "$PY" ]]; then
    log "ERROR: no python3 found in PATH"
    exit 3
fi
log "python: $PY ($($PY --version 2>&1))"

[[ -d scripts ]] || { log "ERROR: scripts/ missing"; exit 4; }

hdr "paleotopo-data-assimilation pipeline started"
log "  no-videos=$NO_VIDEOS  no-pygmt=$NO_PYGMT  no-sensitivity=$NO_SENSITIVITY  no-thickness=$NO_THICKNESS  force=$FORCE"
log "  pygmt cadence=$PYGMT_CADENCE Ma"
T0=$(date +%s)

FORCE_ARG=""
RENDER_FORCE_ARG=""
if [[ $FORCE -eq 1 ]]; then
    hdr "FORCE cleanup"
    log "  removing plate-ID cache"
    rm -fv data/geochem/sample_plate_ids_SW.npy 2>/dev/null
    log "  removing video frame caches (they'd be re-used otherwise)"
    rm -rfv outputs/video_frames_cartopy_SW outputs/video_frames_pygmt_SW 2>/dev/null
    FORCE_ARG="--force"
    RENDER_FORCE_ARG="--force"
fi

pushd scripts > /dev/null

log "  [1/6] assimilating all 109 slices …"
run_tagged "assim"      $PY -u assimilate_scotese.py --all $FORCE_ARG

log "  [2/6] building summary stats …"
run_tagged "stats"      $PY -u build_summary_stats_scotese.py

log "  [3/6] full-sweep diagnostic figures (Figs 6, 7, 9) …"
run_tagged "diag"       $PY -u full_sweep_diagnostics_scotese.py

if [[ $NO_SENSITIVITY -eq 0 ]]; then
    log "  [4/7] four-refinement ablation (Fig 10) …"
    run_tagged "ablation" $PY -u sensitivity_refinements.py
else
    log "  [4/7] skipping sensitivity ablation (--no-sensitivity)"
fi

if [[ $NO_THICKNESS -eq 0 ]]; then
    log "  [5/7] Airy-isostasy crustal-thickness derivation (Fig 11) …"
    run_tagged "thickness" $PY -u derive_crustal_thickness.py --all --figure
else
    log "  [5/7] skipping crustal-thickness derivation (--no-thickness)"
fi

if [[ $NO_VIDEOS -eq 0 ]]; then
    log "  [6/7] cartopy preview videos (5 modes: 4 elevation + crustal thickness) …"
    run_tagged "cartopy" $PY -u render_videos_cartopy_scotese.py $RENDER_FORCE_ARG
    if [[ $NO_THICKNESS -eq 0 ]]; then
        log "         + crustal-thickness MP4 …"
        run_tagged "thickness.video" $PY -u render_video_crustal_thickness.py $RENDER_FORCE_ARG
    fi
    if [[ $NO_PYGMT -eq 0 ]]; then
        log "  [7/7] pyGMT publication videos (Winkel-Tripel, cadence $PYGMT_CADENCE Ma) …"
        run_tagged "pygmt" $PY -u render_videos_pygmt_scotese.py \
            --cadence "$PYGMT_CADENCE" --fps 10 $RENDER_FORCE_ARG
    else
        log "  [7/7] skipping pyGMT videos (--no-pygmt-videos)"
    fi
else
    log "  [6-7/7] skipping ALL videos (--no-videos)"
fi

popd > /dev/null

T1=$(date +%s); DUR=$((T1 - T0))
HOURS=$((DUR / 3600)); MINUTES=$(( (DUR % 3600) / 60 ))
hdr "pipeline complete in ${HOURS}h ${MINUTES}m"
log "Outputs:"
log "  Corrected NetCDFs : data/corrected/"
log "  Crustal-thickness NetCDFs : data/corrected/<age>Ma_crustal_thickness_SW.nc"
log "  Stats + figures   : data/corrected/  +  outputs/"
log "  Paper figures     : Figures/Fig01..Fig11_*  (run ./build_all_figures.sh"
log "                       to assemble them all into Figures/)"
log "  Videos            : outputs/SW_paleotopo_*.mp4"
