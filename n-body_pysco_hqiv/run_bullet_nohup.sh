#!/usr/bin/env bash
# Run Bullet simulation detached and (optionally) with resource limits so it
# cannot starve Cursor / the rest of the system. Use for heavy runs.

set -e
cd "$(dirname "$0")"

OUT="${1:-./bullet_5mpc_512/}"
LOG="${2:-bullet_5mpc_512.log}"
shift 2 2>/dev/null || true
# Only pass non-empty extra args (avoid "unrecognized arguments: " from empty string)
EXTRA_ARGS=()
for a in "$@"; do [ -n "$a" ] && EXTRA_ARGS+=("$a"); done

# Optional: cap memory [GB] and CPU [%] so the run can't kill the machine.
# Example: BULLET_MEM_GB=8 BULLET_CPU_PCT=70 ./run_bullet_nohup.sh
BULLET_MEM_GB="${BULLET_MEM_GB:-}"
BULLET_CPU_PCT="${BULLET_CPU_PCT:-}"

# Use all cores by default (set before python so BLAS/FFT see it)
NTHREADS="${BULLET_NTHREADS:-$(nproc 2>/dev/null || echo 32)}"
export OMP_NUM_THREADS="$NTHREADS"
export OPENBLAS_NUM_THREADS="$NTHREADS"
export MKL_NUM_THREADS="$NTHREADS"
export NUMEXPR_NUM_THREADS="$NTHREADS"

_run() {
  # Optional: limit virtual memory (ulimit -v is in KB); no systemd needed
  if [ -n "$BULLET_MEM_GB" ]; then
    ulimit -v $((BULLET_MEM_GB * 1024 * 1024)) 2>/dev/null || true
  fi
  nice -n 10 nohup python3 -u run_bullet.py \
    --resolution 512 \
    --npart 1e7 \
    --box 5 \
    --output "$OUT" \
    --nthreads "$NTHREADS" \
    "${EXTRA_ARGS[@]}" \
    >> "$LOG" 2>&1 &
  echo $!
}

if command -v systemd-run >/dev/null 2>&1 && { [ -n "$BULLET_MEM_GB" ] || [ -n "$BULLET_CPU_PCT" ]; }; then
  # Run inside a scope with limits (so it can't take down Cursor)
  ARGS=()
  [ -n "$BULLET_MEM_GB" ] && ARGS+=( -p "MemoryMax=${BULLET_MEM_GB}G" )
  [ -n "$BULLET_CPU_PCT" ] && ARGS+=( -p "CPUQuota=${BULLET_CPU_PCT}%" )
  EXTRA=$(printf '%q ' "${EXTRA_ARGS[@]}")
  systemd-run --user --scope "${ARGS[@]}" -- env OMP_NUM_THREADS="$NTHREADS" OPENBLAS_NUM_THREADS="$NTHREADS" MKL_NUM_THREADS="$NTHREADS" bash -c "cd '$PWD' && nohup python3 -u run_bullet.py --resolution 512 --npart 1e7 --box 5 --output '$OUT' --nthreads $NTHREADS $EXTRA >> '$LOG' 2>&1" &
  PID=$!
  disown -h $PID
  echo "Started in resource-limited scope (PID=$PID)"
  echo "  MemoryMax=${BULLET_MEM_GB}G  CPUQuota=${BULLET_CPU_PCT}%"
else
  PID=$(_run "${EXTRA_ARGS[@]}")
  disown -h $PID
  echo "Started detached run: PID=$PID"
  if [ -n "$BULLET_MEM_GB" ]; then
    echo "  (ulimit -v ${BULLET_MEM_GB}G; systemd not used)"
  fi
  if [ -z "$BULLET_MEM_GB" ] && [ -z "$BULLET_CPU_PCT" ]; then
    echo ""
    echo "  WARNING: No resource limits. Heavy runs can freeze or kill Cursor."
    echo "  To cap usage, run instead:"
    echo "    BULLET_MEM_GB=8 BULLET_CPU_PCT=70 $0 $OUT $LOG"
    echo "  Or run the job on another machine / job queue."
  fi
fi

echo "  nthreads: $NTHREADS"
echo "  output: $OUT"
echo "  log:    $LOG"
echo "  tail -f $LOG"
