#!/usr/bin/env bash
# run_magnav_sim.sh — Launch PX4 SITL + MagNav magnetic map bridge
#
# Usage:
#   bash gazebo/plugins/run_magnav_sim.sh           # full pipeline
#   bash gazebo/plugins/run_magnav_sim.sh --no-tl   # anomaly map only
#   bash gazebo/plugins/run_magnav_sim.sh --no-anomaly  # TL only
#
# Prerequisites:
#   1. julia gazebo/plugins/export_ottawa_map.jl   (once)
#   2. PX4-Autopilot checkout in ~/PX4-Autopilot
#   3. pip install numpy scipy pyulog
#
# After the flight:
#   python3 gazebo/scripts/test_tl_on_sim.py <path/to/flight.ulg>

set -e
REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
BRIDGE="$REPO_ROOT/gazebo/plugins/magnetic_map_bridge.py"
MAP="$REPO_ROOT/gazebo/plugins/data/ottawa_map.npz"
PX4_DIR="${PX4_DIR:-$HOME/PX4-Autopilot}"

# Ottawa test area coordinates
export PX4_HOME_LAT="${PX4_HOME_LAT:-44.76}"
export PX4_HOME_LON="${PX4_HOME_LON:-76.15}"   # PX4 uses positive West for lon
export PX4_HOME_ALT="${PX4_HOME_ALT:-120}"

# Add our model to Gazebo resource path
export GZ_SIM_RESOURCE_PATH="$REPO_ROOT/gazebo/models:${GZ_SIM_RESOURCE_PATH:-}"

echo "══════════════════════════════════════════════════"
echo "  MagNav SITL Pipeline"
echo "══════════════════════════════════════════════════"
echo "  Ottawa: lat=$PX4_HOME_LAT  lon=-$PX4_HOME_LON"
echo "  Model path: $REPO_ROOT/gazebo/models"

# ── Step 1: Ensure map is exported ────────────────────────────────────────────
if [ ! -f "$MAP" ]; then
    echo ""
    echo "▶ Ottawa map NPZ not found. Exporting from MagNav.jl..."
    julia "$REPO_ROOT/gazebo/plugins/export_ottawa_map.jl"
fi
echo "  Map: $MAP  ✓"

# ── Step 2: Start bridge in background ────────────────────────────────────────
BRIDGE_ARGS="$@"
echo ""
echo "▶ Starting magnetic map bridge..."
python3 "$BRIDGE" $BRIDGE_ARGS &
BRIDGE_PID=$!
echo "  Bridge PID: $BRIDGE_PID"

# Cleanup bridge on exit
cleanup() {
    echo ""
    echo "▶ Stopping bridge (PID=$BRIDGE_PID)..."
    kill "$BRIDGE_PID" 2>/dev/null || true
}
trap cleanup EXIT

# ── Step 3: Launch PX4 SITL ───────────────────────────────────────────────────
echo ""
echo "▶ Launching PX4 SITL (Gazebo Harmonic, x500_magnav)..."
echo "  PX4_DIR: $PX4_DIR"
echo ""
echo "  [!] When simulation starts, load the box mission:"
echo "      gazebo/scripts/box_mission.plan"
echo ""

cd "$PX4_DIR"
make px4_sitl gz_x500_magnav

echo ""
echo "══════════════════════════════════════════════════"
echo "  Simulation complete."
echo ""
echo "  Next steps:"
echo "  1. Find the .ulg log:"
echo "     ls -lt $PX4_DIR/build/px4_sitl_default/logs/*.ulg | head -1"
echo "  2. Run TL validation:"
echo "     python3 $REPO_ROOT/gazebo/scripts/test_tl_on_sim.py <flight.ulg>"
echo "══════════════════════════════════════════════════"
