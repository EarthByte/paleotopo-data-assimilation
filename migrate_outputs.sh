#!/usr/bin/env bash
# =============================================================================
# migrate_outputs.sh — one-off cleanup of outputs[_Scotese]/ directory layout
# =============================================================================
#
# Moves existing PNG/PDF/CSV files to outputs[_Scotese]/figures/
# and MP4 videos to outputs[_Scotese]/videos/, matching the new layout
# the figure/video scripts now write to.
#
# Frame caches (video_frames_*) stay at the top level — they're
# temporary per-frame PNG dumps that get wiped between runs.
#
# Safe to re-run.  Files already in the right place stay there.
#
# USAGE
#   cd <project>/
#   bash migrate_outputs.sh
# =============================================================================
set -u

migrate_one() {
    local DIR="$1"           # outputs/ or outputs_Scotese/
    [[ ! -d "$DIR" ]] && { echo "skip: $DIR does not exist"; return; }
    echo ""
    echo "=== $DIR ==="
    mkdir -p "$DIR/figures" "$DIR/videos"

    # Move static outputs (PNG/PDF/CSV) into figures/
    local moved_fig=0
    for f in "$DIR"/*.png "$DIR"/*.pdf "$DIR"/*.csv; do
        [[ ! -e "$f" ]] && continue
        # Don't move files already under figures/ or videos/
        local base
        base="$(basename "$f")"
        if [[ -e "$DIR/figures/$base" ]]; then
            # Overwrite-with-warn: prefer the newer file
            if [[ "$f" -nt "$DIR/figures/$base" ]]; then
                mv -f "$f" "$DIR/figures/$base"
                echo "  moved (replaced older): $base → figures/"
                moved_fig=$((moved_fig + 1))
            else
                echo "  skipped (newer in figures/): $base"
                rm -f "$f"  # remove the older duplicate at top level
            fi
        else
            mv "$f" "$DIR/figures/"
            echo "  moved: $base → figures/"
            moved_fig=$((moved_fig + 1))
        fi
    done

    # Move videos (MP4) into videos/
    local moved_vid=0
    for f in "$DIR"/*.mp4; do
        [[ ! -e "$f" ]] && continue
        local base
        base="$(basename "$f")"
        if [[ -e "$DIR/videos/$base" ]]; then
            if [[ "$f" -nt "$DIR/videos/$base" ]]; then
                mv -f "$f" "$DIR/videos/$base"
                echo "  moved (replaced older): $base → videos/"
                moved_vid=$((moved_vid + 1))
            else
                echo "  skipped (newer in videos/): $base"
                rm -f "$f"
            fi
        else
            mv "$f" "$DIR/videos/"
            echo "  moved: $base → videos/"
            moved_vid=$((moved_vid + 1))
        fi
    done

    echo "  $moved_fig figure(s) and $moved_vid video(s) moved"
    echo "  frame caches left in place: $(ls -d "$DIR"/video_frames_* 2>/dev/null | wc -l | tr -d ' ') directories"
}

cd "$(dirname "$0")"
migrate_one "outputs"
migrate_one "outputs_Scotese"

echo ""
echo "Done.  New layout:"
echo "  outputs/figures/         outputs_Scotese/figures/"
echo "  outputs/videos/          outputs_Scotese/videos/"
echo "  outputs/video_frames_*/  outputs_Scotese/video_frames_*/   (caches, unchanged)"
