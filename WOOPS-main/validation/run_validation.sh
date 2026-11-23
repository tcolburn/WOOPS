#!/usr/bin/env bash
# Execute both WOOPS programs on the committed real wildtype excerpt, then compare.
set -euo pipefail
VAL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${VAL_DIR}/.." && pwd)"
BUILD=1
if [[ "${1:-}" == '--no-build' ]]; then BUILD=0; shift; fi
if [[ "$#" -ne 0 ]]; then echo 'Usage: validation/run_validation.sh [--no-build]' >&2; exit 2; fi
for branch in nanodomain kirkwood; do
  if [[ ! -d "$REPO_DIR/src/$branch" ]]; then
    echo "ERROR: missing $REPO_DIR/src/$branch; put validation/ in repo root." >&2
    exit 2
  fi
done
if [[ "$BUILD" == 1 ]]; then
  make -C "$REPO_DIR/src/nanodomain"
  make -C "$REPO_DIR/src/kirkwood"
fi
for branch in nanodomain kirkwood; do
  exe="$REPO_DIR/src/$branch/map.e"
  [[ -x "$exe" ]] || { echo "ERROR: missing executable $exe" >&2; exit 2; }
  work="$VAL_DIR/run/$branch"
  mkdir -p "$work"
  # Remove stale results so a failed run cannot accidentally pass validation.
  if [[ "$branch" == 'nanodomain' ]]; then
    rm -f "$work/p1p2.dat" "$work/restart"
  else
    rm -f "$work"/kirkwood_*.dat "$work/restart"
  fi
  echo "Running WOOPS $branch on real wildtype excerpt"
  (cd "$work" && OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}" "$exe" < WOOPS_input > run.log 2>&1) || {
     echo "ERROR: $branch exited unsuccessfully. See $work/run.log" >&2
     exit 1
  }
  echo "Completed $branch; log: $work/run.log"
done
python "$VAL_DIR/validate.py"
