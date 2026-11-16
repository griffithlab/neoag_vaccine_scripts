#!/usr/bin/env bash
# Run the three genomic review QC scripts with one command:
#   get_neoantigen_qc.py, get_FDA_thresholds.py, hla_comparison.py
#
# Usage (inside the griffithlab/neoang_scripts container):
#   bash /opt/scripts/run_manual_review.sh -w $WORKING_BASE -y $HOME/yamls/${GCS_CASE_NAME}_immuno_cloud-WDL.yaml
#
# Options:
#   -w  WORKING_BASE (gcp_immuno folder)                       [required]
#   -y  pipeline YAML file (needed by get_neoantigen_qc.py)    [required]
#   -f  final results folder name                              [default: final_results]
#   -o  output folder                                          [default: WORKING_BASE/../manual_review]
#
# Each script runs even if an earlier one fails; the exit code is non-zero if any failed.

set -u

usage() {
    sed -n '2,14p' "$0" | sed 's/^# \{0,1\}//'
    exit 1
}

WB=""
YAML=""
FIN="final_results"
OUT=""

while getopts "w:y:f:o:h" opt; do
    case "$opt" in
        w) WB="$OPTARG" ;;
        y) YAML="$OPTARG" ;;
        f) FIN="$OPTARG" ;;
        o) OUT="$OPTARG" ;;
        *) usage ;;
    esac
done

[ -z "$WB" ]   && { echo "ERROR: -w WORKING_BASE is required" >&2; usage; }
[ -z "$YAML" ] && { echo "ERROR: -y YAML file is required" >&2; usage; }
[ -d "$WB" ]   || { echo "ERROR: WORKING_BASE not found: $WB" >&2; exit 1; }
[ -f "$YAML" ] || { echo "ERROR: YAML file not found: $YAML" >&2; exit 1; }

OUT="${OUT:-$WB/../manual_review}"
mkdir -p "$OUT"

# Find the python scripts next to this one (works in /opt/scripts or a local checkout)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="${PYTHON:-python3}"

FAILED=()

run_step() {
    local name="$1"; shift
    echo
    echo "################################################################################"
    echo "# $name"
    echo "################################################################################"
    if ! "$PY" "$SCRIPT_DIR/$name" "$@"; then
        echo "!! $name FAILED" >&2
        FAILED+=("$name")
    fi
}

run_step get_neoantigen_qc.py -WB "$WB" -f "$FIN" --yaml "$YAML" --o "$OUT"
run_step get_FDA_thresholds.py -WB "$WB" -f "$FIN" --o "$OUT"
run_step hla_comparison.py     -WB "$WB" -f "$FIN" --o "$OUT"

echo
echo "Outputs written to: $(cd "$OUT" && pwd)"
if [ ${#FAILED[@]} -gt 0 ]; then
    echo "Failed scripts: ${FAILED[*]}" >&2
    exit 1
fi
echo "All three scripts completed."
