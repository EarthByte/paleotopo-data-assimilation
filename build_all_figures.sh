#!/usr/bin/env bash
# =============================================================================
# build_all_figures.sh  —  rebuild every paper-figure script (Scotese)
# =============================================================================
#
# Runs every figure-producing script in scripts/ in dependency
# order.  Use this after a fresh Scotese assimilation sweep, or any time
# the Scotese knobs change and you want a clean set of figures.
#
# Prereq: data/corrected/<age>Ma_corrected.nc must exist for the
# ages each figure script needs.  Build them with the Scotese assimilation
# runner first (the Scotese equivalent of run_pipeline.sh).
#
# USAGE
#   cd <repo root>
#   ./build_all_figures.sh                 # build everything
#   ./build_all_figures.sh --no-pygmt      # skip the pyGMT comparison
#                                                 # figure (use cartopy only)
#   ./build_all_figures.sh -h | --help     # show this help
#
# OUTPUT
#   Figures/Fig01..Fig10*.png        — paper-numbered figures
#                                            (single source of truth for the
#                                             Earth-Science Reviews paper)
#   outputs/SW_*.png               — intermediate diagnostics +
#                                            previews used as inputs to the
#                                            paper-numbered copies above
#   data/corrected/*.csv, *.md     — refreshed summary stats
#
# =============================================================================
set -u
set -o pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

WITH_PYGMT=1
while [[ $# -gt 0 ]]; do
    case "$1" in
        --no-pygmt)  WITH_PYGMT=0 ;;
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
# Fail-stop wrapper (see run_pipeline.sh for the rationale).
run_or_die() {
    local tag="$1"
    run_tagged "$@"
    local rc=$?
    if [[ $rc -ne 0 ]]; then
        hdr "STAGE FAILED — [$tag] exited with code $rc"
        log "Aborting figure build.  Fix the error above and re-run."
        exit "$rc"
    fi
}

PY=$(command -v python3 || command -v python)
[[ -n "$PY" ]] || { log "ERROR: no python3 found in PATH"; exit 3; }
[[ -d scripts ]] || { log "ERROR: scripts/ missing"; exit 4; }
log "python: $PY ($($PY --version 2>&1))"

hdr "build_all_figures.sh started"
log "  with-pygmt-comparison=$WITH_PYGMT"
T0=$(date +%s)

PROJ_ROOT="$HERE"
PFIG="$PROJ_ROOT/Figures"
OUTS="$PROJ_ROOT/outputs"
mkdir -p "$PFIG"

pushd scripts > /dev/null

# Step 0 — refresh summary CSVs (Scotese versions have the _scotese suffix).
log "  [0/10] refresh summary stats CSVs + markdown dashboard …"
run_or_die "stats"    $PY -u build_summary_stats_scotese.py

# Step 1 — methodology flowchart.  Script writes Fig01 straight into
# Figures/ (plus a mirror copy in outputs/).
log "  [1/10] methodology flowchart (Fig01) …"
run_or_die "flow"     $PY -u draw_methodology_flowchart.py

# Step 2 — sample-distribution map (Fig02 → Figures/ directly).
log "  [2/10] sample-distribution figure (Fig02) …"
run_or_die "samples"  $PY -u make_sample_distribution_figure.py

# Step 3 — single-slice pipeline illustration (Fig03 → Figures/ directly).
log "  [3/10] pipeline-illustration figure (Fig03) …"
run_or_die "illustr"  $PY -u make_pipeline_illustration_figure.py

# Step 4 — full-sweep diagnostics.  Script writes SW_*.png to
# outputs/; copy to Figures/ with Fig06/07/09 names.
log "  [4/10] full-sweep diagnostic figures (Fig06/07/09) …"
run_or_die "diag"     $PY -u full_sweep_diagnostics_scotese.py
cp "$OUTS/SW_full_sweep_diagnostics.png"   "$PFIG/Fig06_full_sweep_diagnostics.png"
cp "$OUTS/SW_hypsometry_selected_ages.png" "$PFIG/Fig07_hypsometric_curves.png"
cp "$OUTS/SW_metrics_by_era.png"           "$PFIG/Fig09_metrics_by_era.png"
log "    copied Fig06/07/09 → $PFIG/"

# Step 5 — Fig05a/b/c comparison figures.
#
# Two renderers are available:
#   - make_comparison_figures.py         (pyGMT, Winkel-Tripel — publication)
#   - make_comparison_figures_cartopy.py (cartopy, Robinson    — fallback)
#
# Both produce the same Fig05a/b/c.png filenames in Figures/, so we
# only run one.  Default: pyGMT.  --no-pygmt switches to the cartopy
# fallback (useful if GMT isn't installed on this machine).
if [[ $WITH_PYGMT -eq 1 ]]; then
    log "  [5/10] comparison figures (pyGMT, Winkel-Tripel) (Fig05a/b/c) …"
    run_or_die "cmp"  $PY -u make_comparison_figures.py
    for pair in 500-400 300-200 100-50; do
        label="a"; case $pair in 300-200) label="b" ;; 100-50) label="c" ;; esac
        cp "$OUTS/SW_comparison_${pair}Ma.png" \
           "$PFIG/Fig05${label}_SW_comparison_${pair}Ma.png"
        # PDF is optional — only copy if pyGMT produced one
        if [[ -f "$OUTS/SW_comparison_${pair}Ma.pdf" ]]; then
            cp "$OUTS/SW_comparison_${pair}Ma.pdf" \
               "$PFIG/Fig05${label}_SW_comparison_${pair}Ma.pdf"
        fi
    done
    log "    copied Fig05a/b/c (pyGMT) → $PFIG/"
else
    log "  [5/10] comparison figures (cartopy fallback) (Fig05a/b/c) …"
    run_or_die "cmp-cp"   $PY -u make_comparison_figures_cartopy.py
fi

# Step 6 — supermountain epochs (Fig08 → Figures/ directly).
log "  [6/10] supermountains figure (Fig08) …"
run_or_die "supermt"  $PY -u make_supermountains_figure.py

# Step 7 — plate-ID validation (Fig04 → Figures/ directly).
log "  [7/10] plate-ID validation figure (Fig04) …"
run_or_die "plateid"  $PY -u make_plate_id_validation_figure.py

# Step 8 — correction-distance diagnostic (not a paper figure; lands in
# outputs/).
log "  [8/10] correction-distance diagnostic …"
run_or_die "distdiag" $PY -u diagnose_correction_distances.py

# Step 9 — knob-sweep sensitivity figure (Fig10 → Figures/ directly).
log "  [9/10] sensitivity-refinement figure (Fig10) …"
run_or_die "sens"     $PY -u sensitivity_refinements.py

popd > /dev/null

T1=$(date +%s); DUR=$((T1 - T0))
MINUTES=$((DUR / 60))
SECONDS=$((DUR % 60))
hdr "all Scotese figures rebuilt in ${MINUTES}m ${SECONDS}s"

# Show the user the final paper-figures inventory — that's the single
# source of truth for what the manuscript embeds.
log "Figures/ inventory (sorted):"
ls -1 Figures/ 2>/dev/null | sort | sed 's|^|    |'
