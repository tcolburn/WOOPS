#!/usr/bin/env bash
# Run the *actual* WOOPS executables on the included wildtype trajectory.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
BUILD=1
COMPARE=1
for arg in "$@"; do
  case "$arg" in
    --no-build) BUILD=0 ;;
    --no-compare) COMPARE=0 ;;
    -h|--help)
      echo 'Usage: bash examples/wildtype_plastocyanin/run_example.sh [--no-build] [--no-compare]'
      exit 0 ;;
    *) echo "Unknown option: $arg" >&2; exit 2 ;;
  esac
done
if [ ! -f "$REPO/src/nanodomain/makefile" ] || [ ! -f "$REPO/src/kirkwood/makefile" ]; then
  echo "Missing WOOPS source branches under $REPO/src/" >&2
  echo 'Place examples/wildtype_plastocyanin/ inside the WOOPS repository.' >&2
  exit 2
fi
: "${OMP_NUM_THREADS:=2}"
export OMP_NUM_THREADS
if [ "$BUILD" -eq 1 ]; then
  echo '[1/3] Building WOOPS (GNU Fortran/OpenMP) ...'
  make -C "$REPO/src/nanodomain"
  make -C "$REPO/src/kirkwood"
else
  echo '[1/3] Skipping compilation (--no-build)'
fi
for branch in nanodomain kirkwood; do
  exe="$REPO/src/$branch/map.e"
  if [ ! -x "$exe" ]; then
    echo "Executable missing: $exe" >&2
    exit 2
  fi
  rundir="$HERE/runs/$branch"
  # Remove *only* results of the prior example run, never the saved examples.
  if [ "$branch" = 'nanodomain' ]; then
    rm -f "$rundir/p1p2.dat"
  else
    rm -f "$rundir"/kirkwood_*.dat
  fi
  rm -f "$rundir/restart" "$rundir/restart."* "$rundir/run.log"
  echo "[2/3] Running $branch on real wildtype DCD (10 frames) ..."
  ( cd "$rundir" && "$exe" < WOOPS_input > run.log 2>&1 ) || {
    echo "$branch failed; see $rundir/run.log" >&2
    tail -n 35 "$rundir/run.log" >&2 || true
    exit 1
  }
  if ! grep -q 'COMPLETED SUCCESSFULLY' "$rundir/run.log"; then
    echo "$branch did not report completion; see $rundir/run.log" >&2
    exit 1
  fi
done
if [ "$COMPARE" -eq 1 ]; then
  echo '[3/3] Comparing newly generated data to the original WOOPS results ...'
  python3 "$HERE/scripts/compare_results.py"
else
  echo '[3/3] Finished (--no-compare). Outputs are in runs/{nanodomain,kirkwood}/.'
fi
