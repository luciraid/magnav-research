#!/usr/bin/env bash
# run_test.sh — End-to-end TL validation on a Gazebo SITL flight.
#
# USAGE
#   ./gazebo/run_test.sh <coef_file> [--ulg <path/to/flight.ulg>]
#
# MODES
#   Interactive (no --ulg):
#     Launches PX4 SITL at Ottawa, waits for you to fly a mission in QGC,
#     prompts you to land, then auto-detects the newest .ulg and runs the test.
#
#   Direct (--ulg provided):
#     Skips PX4 launch entirely.  Runs the test immediately on the given log.
#     Use this to re-run analysis on a log you already have.
#
# EXAMPLES
#   ./gazebo/run_test.sh gazebo/test_data/zero_coef.txt
#   ./gazebo/run_test.sh gazebo/test_data/my_coef.txt
#   ./gazebo/run_test.sh gazebo/test_data/my_coef.txt --ulg /path/to/18_01_09.ulg
#
# OUTPUT
#   gazebo/results/<coef_stem>_<timestamp>/
#     metrics.txt        — full stdout from test_tl_on_sim.py
#     <log>_tl_validation.png — 4-panel summary plot
#
# AFTER BOTH TESTS
#   python3 gazebo/scripts/generate_comparison.py \
#       gazebo/results/zero_<ts>  gazebo/results/my_<ts>
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

# ── Repo root detection ───────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

# ── PX4 log dir ───────────────────────────────────────────────────────────────
PX4_LOG_DIR="${HOME}/PX4-Autopilot/build/px4_sitl_default/rootfs/log"

# ── Defaults ─────────────────────────────────────────────────────────────────
COEF_FILE=""
ULG_FILE=""
LAUNCH_PX4=true

# ── Argument parsing ─────────────────────────────────────────────────────────
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <coef_file> [--ulg <flight.ulg>]"
    echo ""
    echo "  coef_file  : TL coefficient file (e.g. gazebo/test_data/zero_coef.txt)"
    echo "               If just a filename, looks in gazebo/test_data/ automatically."
    echo "  --ulg PATH : Skip PX4 launch, run directly on this .ulg file."
    exit 1
fi

COEF_ARG="$1"
shift

while [[ $# -gt 0 ]]; do
    case "$1" in
        --ulg)
            ULG_FILE="$2"
            LAUNCH_PX4=false
            shift 2
            ;;
        *)
            echo "Unknown argument: $1"
            exit 1
            ;;
    esac
done

# ── Resolve coef file ─────────────────────────────────────────────────────────
if [[ -f "$COEF_ARG" ]]; then
    COEF_FILE="$(realpath "$COEF_ARG")"
elif [[ -f "${REPO_ROOT}/gazebo/test_data/${COEF_ARG}" ]]; then
    COEF_FILE="$(realpath "${REPO_ROOT}/gazebo/test_data/${COEF_ARG}")"
elif [[ -f "${REPO_ROOT}/gazebo/test_data/$(basename "$COEF_ARG")" ]]; then
    COEF_FILE="$(realpath "${REPO_ROOT}/gazebo/test_data/$(basename "$COEF_ARG")")"
else
    echo "ERROR: Cannot find coefficient file: $COEF_ARG"
    echo "  Looked in: current dir, ${REPO_ROOT}/gazebo/test_data/"
    exit 1
fi

COEF_STEM="$(basename "${COEF_FILE%.txt}")"

# ── Output directory ──────────────────────────────────────────────────────────
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
OUT_DIR="${REPO_ROOT}/gazebo/results/${COEF_STEM}_${TIMESTAMP}"
mkdir -p "$OUT_DIR"

echo "========================================================"
echo "  MagNav TL Validation — run_test.sh"
echo "========================================================"
echo ""
echo "  Coefficients : $COEF_FILE"
echo "  Output dir   : $OUT_DIR"
echo ""

# ── Optional: Launch PX4 SITL ─────────────────────────────────────────────────
if $LAUNCH_PX4; then
    echo "── Step 1: Launch PX4 SITL ─────────────────────────────"
    echo ""
    echo "  Home position: Ottawa area"
    echo "    LAT=44.76  LON=-76.15  ALT=397"
    echo ""

    # Snapshot of current newest log (to detect new log after flight)
    SNAPSHOT_TIME=$(date +%s)

    echo "  Launching PX4 SITL in background..."
    echo "  (Set GZ_SIM_SYSTEM_PLUGIN_PATH if using the C++ MagneticMapPlugin)"
    echo ""

    # Export plugin path if built
    PLUGIN_SO="${REPO_ROOT}/build/plugin/libMagneticMapPlugin.so"
    if [[ -f "$PLUGIN_SO" ]]; then
        export GZ_SIM_SYSTEM_PLUGIN_PATH="$(dirname "$PLUGIN_SO")"
        echo "  Plugin found — GZ_SIM_SYSTEM_PLUGIN_PATH=$GZ_SIM_SYSTEM_PLUGIN_PATH"
    fi

    (
        cd "${HOME}/PX4-Autopilot"
        PX4_HOME_LAT=44.76 PX4_HOME_LON=-76.15 PX4_HOME_ALT=397 \
            make px4_sitl gz_x500 2>&1 | tee "${OUT_DIR}/px4_sitl.log"
    ) &
    PX4_PID=$!

    echo "  PX4 PID: $PX4_PID"
    echo ""
    echo "  ▶ Open QGroundControl and fly your mission."
    echo "  ▶ After landing and disarming, press ENTER here to continue."
    echo ""
    read -r -p "  Press ENTER after landing (or Ctrl+C to abort): "
    echo ""

    echo "── Step 2: Stop PX4 ────────────────────────────────────"
    kill "$PX4_PID" 2>/dev/null || true
    sleep 2
    echo ""

    echo "── Step 3: Auto-detect newest .ulg ─────────────────────"
    ULG_FILE=$(find "$PX4_LOG_DIR" -name "*.ulg" \
        -newer <(touch -t "$(date -d @$SNAPSHOT_TIME +%Y%m%d%H%M.%S)" /tmp/_ts_ref && echo /tmp/_ts_ref | xargs) \
        -printf '%T@ %p\n' 2>/dev/null | sort -rn | head -1 | awk '{print $2}')

    if [[ -z "$ULG_FILE" ]]; then
        # Fallback: just newest in log dir
        ULG_FILE=$(find "$PX4_LOG_DIR" -name "*.ulg" -printf '%T@ %p\n' 2>/dev/null \
            | sort -rn | head -1 | awk '{print $2}')
    fi

    if [[ -z "$ULG_FILE" ]]; then
        echo "ERROR: No .ulg file found in $PX4_LOG_DIR"
        exit 1
    fi
    echo "  Newest log: $ULG_FILE"
    echo ""
else
    # ── Direct mode: validate the provided .ulg ───────────────────────────────
    if [[ ! -f "$ULG_FILE" ]]; then
        echo "ERROR: .ulg file not found: $ULG_FILE"
        exit 1
    fi
    ULG_FILE="$(realpath "$ULG_FILE")"
    echo "── Using existing flight log ────────────────────────────"
    echo "  $ULG_FILE"
    echo "  Size: $(du -h "$ULG_FILE" | cut -f1)"
    echo ""
fi

# Save the log path for later (generate_comparison.py reads this)
echo "$ULG_FILE" > "${OUT_DIR}/ulg_path.txt"
echo "$COEF_FILE" > "${OUT_DIR}/coef_path.txt"

# ── Run test_tl_on_sim.py ─────────────────────────────────────────────────────
echo "── Running TL validation ────────────────────────────────"
echo "  python3 gazebo/scripts/test_tl_on_sim.py"
echo "          $ULG_FILE"
echo "          --coef $COEF_FILE"
echo "          --output $OUT_DIR"
echo ""

python3 "${REPO_ROOT}/gazebo/scripts/test_tl_on_sim.py" \
    "$ULG_FILE" \
    --coef "$COEF_FILE" \
    --output "$OUT_DIR" \
    | tee "${OUT_DIR}/metrics.txt"

PLOT=$(find "$OUT_DIR" -name "*.png" | head -1)

echo ""
echo "========================================================"
echo "  ✓ Test complete"
echo "  Results : $OUT_DIR"
echo "  Metrics : ${OUT_DIR}/metrics.txt"
if [[ -n "$PLOT" ]]; then
    echo "  Plot    : $PLOT"
fi
echo ""
echo "  To generate a comparison report after running both tests:"
echo "    python3 gazebo/scripts/generate_comparison.py \\"
echo "        <zero_results_dir>  <my_results_dir>"
echo "========================================================"
