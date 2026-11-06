#!/usr/bin/env bash
# =============================================================================
# build_all_figures.sh  —  rebuild the figures and tables of the paper (Scotese)
# =============================================================================
#
# Runs every figure and table script in scripts/ in dependency order and
# gives each figure its number in the Earth-Science Reviews paper.
#
# Prereq: the 109 corrected grids data/corrected/<age>Ma_corrected_SW.nc
# (from the Zenodo archive, or rebuilt with ./run_pipeline.sh) and the
# inputs listed in data/README.md. The regional comparisons read the
# independent paleoelevation records shipped in data/independent_records/.
#
# USAGE
#   cd <repo root>
#   ./build_all_figures.sh                 # build everything
#   ./build_all_figures.sh --no-pygmt      # skip the pyGMT maps (Figs. 5, 8-12);
#                                          # Tables S4-S8 need them too
#   ./build_all_figures.sh -h | --help     # show this help
#
# OUTPUT
#   Figures/FigNN_<name>.{png,pdf,svg}  the paper-numbered figures:
#       Fig01 sample distribution          Fig09 southern Tibet
#       Fig02 workflow diagram             Fig10 N. American Cordillera
#       Fig03 plate-ID validation          Fig11 Central Pangean Mountains
#       Fig04 workflow at 100 Ma           Fig12 eastern Asia
#       Fig05 input / corrected / dz       Fig13 metrics by era
#       Fig06 full-sweep diagnostics       Fig14 area above 2 and 3 km
#       Fig07 hypsometric curves           Fig15 comparison with the Geology record
#       Fig08 central Andes                Fig16 derived crustal thickness
#       FigS4 sensitivity to the refinements
#       FigS5, FigS6 kernel-length sweep and validation by province (built
#                    only when outputs/kernel_length_sweep/ holds the held-out
#                    cross-validation tables)
#   Figures/tables/TableS4..S8_<region>_records.csv   the independent records
#                    against the input and corrected grids, one table per
#                    region (central Andes, southern Tibet, N. American
#                    Cordillera, Central Pangean Mountains, eastern Asia)
#   outputs/figures/   every figure under its descriptive name
#   outputs/tables/    every table (CSV, and Markdown where written)
#   outputs/esr_revision/   support_km.nc, declustered.csv, per_age.csv
#   data/corrected/*.csv, *.md   refreshed summary statistics
#
# A Figures/Fig* file that this build does not produce (an older numbering)
# is moved to outputs/superseded_paper_figures/, never deleted.
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
# Fail-stop wrapper: a figure built on a failed step is worse than no figure.
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
log "  with-pygmt=$WITH_PYGMT"
T0=$(date +%s)

PFIG="$HERE/Figures"
OUTD="${PALEO_OUTPUT_DIR:-$HERE/outputs}"
OUTS="${PALEO_FIGURES_DIR:-$OUTD/figures}"      # paths_scotese.FIGURES_DIR
TABS="$OUTD/tables"
mkdir -p "$PFIG/tables"

pushd scripts > /dev/null

log "  [1/11] summary statistics (per_slice_stats_SW.csv) …"
run_or_die "stats"    $PY -u build_summary_stats_scotese.py

log "  [2/11] support layer, declustered observations, per-age areas (reproduces the grids to 1 m) …"
run_or_die "diag"     $PY -u esr_revision_diagnostics.py

log "  [3/11] sample distribution (Fig. 1) …"
run_or_die "samples"  $PY -u make_sample_distribution_figure.py

log "  [4/11] workflow diagram (Fig. 2) …"
run_or_die "flow"     $PY -u esr_revision_flowchart.py

log "  [5/11] plate-ID validation, area above thresholds, Geology record (Figs. 3, 14, 15) …"
run_or_die "esrfig"   $PY -u esr_revision_figures.py fig3 fig9 fig10

log "  [6/11] the workflow at 100 Ma (Fig. 4) …"
run_or_die "illustr"  $PY -u make_pipeline_illustration_figure.py

log "  [7/11] full-sweep diagnostics, hypsometry, metrics by era (Figs. 6, 7, 13) …"
run_or_die "sweep"    $PY -u full_sweep_diagnostics_scotese.py

log "  [8/11] derived crustal thickness (Fig. 16) …"
run_or_die "thick"    $PY -u derive_crustal_thickness.py --figure

log "  [9/11] sensitivity to the four refinements (Fig. S4) …"
run_or_die "sens"     $PY -u sensitivity_refinements.py

if [[ -f "$OUTD/kernel_length_sweep/reach_age_summary.csv" && -f "$OUTD/kernel_length_sweep/reach_heldout_samples.csv.gz" ]]; then
    log "  [10/11] kernel-length sweep and validation by province (Figs. S5, S6; Tables S2, S3) …"
    run_or_die "sweepS5" $PY -u esr_revision_figures.py figS5 tableS2
else
    log "  [10/11] skipping Figs. S5, S6 and Tables S2, S3: no held-out tables in $OUTD/kernel_length_sweep/"
fi

if [[ $WITH_PYGMT -eq 1 ]]; then
    log "  [11/11] maps with the support layer (Fig. 5) and the regional comparisons (Figs. 8-12, Tables S4-S8) …"
    run_or_die "maps"    $PY -u esr_revision_maps.py fig5
    run_or_die "regions" $PY -u esr_revision_regions.py all
else
    log "  [11/11] skipping the pyGMT maps (--no-pygmt): Figs. 5, 8-12 and Tables S4-S8 are not rebuilt"
fi

popd > /dev/null

# ---------------------------------------------------------------------------
# Paper numbering: <number> <published name> <file stem in outputs/figures>
# ---------------------------------------------------------------------------
MANIFEST="
Fig01 sample_distribution sample_distribution
Fig02 methodology_flowchart methodology_flowchart
Fig03 plate_id_validation plate_id_validation
Fig04 pipeline_illustration_100Ma pipeline_illustration_100Ma
Fig05 SW_comparison SW_comparison
Fig06 full_sweep_diagnostics SW_full_sweep_diagnostics
Fig07 hypsometric_curves SW_hypsometry_selected_ages
Fig08 region_central_andes region_andes
Fig09 region_southern_tibet region_tibet
Fig10 region_cordillera region_cordillera
Fig11 region_central_pangean_mountains region_cpm
Fig12 region_eastern_asia region_east_asia
Fig13 metrics_by_era SW_metrics_by_era
Fig14 area_above_thresholds area_above_thresholds
Fig15 geology_record_comparison geology_record_comparison
Fig16 derived_crustal_thickness derived_crustal_thickness
FigS4 sensitivity_refinements sensitivity_refinements
FigS5 kernel_length_sweep kernel_length_sweep
FigS6 validation_by_province validation_by_province
"
TABLES="
TableS4 central_andes region_andes_records
TableS5 southern_tibet region_tibet_records
TableS6 cordillera region_cordillera_records
TableS7 central_pangean_mountains region_cpm_records
TableS8 eastern_asia region_east_asia_records
TableS2 kernel_length_sweep kernel_length_sweep_table
TableS3 validation_by_province validation_by_province_table
"

hdr "numbering the figures and tables into Figures/"
expected=" "
n_fig=0
while read -r label name stem; do
    [[ -z "$label" ]] && continue
    for ext in png pdf svg; do
        expected="$expected${label}_${name}.${ext} "
        if [[ -f "$OUTS/$stem.$ext" ]]; then
            cp -p "$OUTS/$stem.$ext" "$PFIG/${label}_${name}.${ext}"
            n_fig=$((n_fig + 1))
        fi
    done
done <<< "$MANIFEST"
n_tab=0
while read -r label name stem; do
    [[ -z "$label" ]] && continue
    for ext in csv md; do
        if [[ -f "$TABS/$stem.$ext" ]]; then
            cp -p "$TABS/$stem.$ext" "$PFIG/tables/${label}_${name}.${ext}"
            n_tab=$((n_tab + 1))
        fi
    done
done <<< "$TABLES"
log "  $n_fig figure files and $n_tab table files written to $PFIG/"

# Figures/Fig* files of an older numbering go to outputs/superseded_paper_figures/, never deleted.
shopt -s nullglob
moved=0
for f in "$PFIG"/Fig*; do
    b=$(basename "$f")
    if [[ "$expected" != *" $b "* ]]; then
        mkdir -p "$OUTD/superseded_paper_figures"
        mv "$f" "$OUTD/superseded_paper_figures/$b"
        moved=$((moved + 1))
    fi
done
shopt -u nullglob
[[ $moved -gt 0 ]] && log "  moved $moved files of an older numbering to $OUTD/superseded_paper_figures/"

T1=$(date +%s); DUR=$((T1 - T0))
MINUTES=$((DUR / 60))
SECONDS=$((DUR % 60))
hdr "figures and tables rebuilt in ${MINUTES}m ${SECONDS}s"

log "Figures/ inventory (sorted):"
ls -1 Figures/ Figures/tables/ 2>/dev/null | sed 's|^|    |'
