# Gazebo MagNav Magnetic Map Plugin

Simulates realistic magnetic anomaly navigation in Gazebo Harmonic + PX4 SITL using the Ottawa area survey map.

## Architecture

```
PX4 SITL (Gazebo Harmonic)
 │
 ├─ /world/.../navsat_sensor/navsat        ─► magnetic_map_bridge.py
 └─ /world/.../magnetometer_sensor/magnetometer ─►     │
                                                         │ 1. Ottawa anomaly (bilinear interp)
                                                         │ 2. 9-term TL interference
                                                         │ 3. Gaussian noise
                                                         ▼
                                              /magnav/magnetometer       ← modified reading
                                              /magnav/mag_anomaly_nT     ← anomaly scalar
                                              /magnav/tl_interference_nT ← TL scalar
```

## Files

| File | Description |
|------|-------------|
| `export_ottawa_map.jl` | Julia: dumps Ottawa map (Eastern_395.h5) → `data/ottawa_map.npz` |
| `magnetic_map_bridge.py` | Python: gz-transport bridge — injects anomaly + TL |
| `data/ottawa_map.npz` | Pre-exported NPZ (created by export step, not in git) |
| `run_magnav_sim.sh` | One-shot launcher: bridge + PX4 SITL |
| `../models/x500_magnav/` | Gazebo model (x500 + zero-noise magnetometer) |

## Signal Model

```
Bt_total = Bt_earth(WMM)  +  anomaly(lat,lon)  +  tl_interference(Bx,By,Bz,c)  +  noise
```

**Ottawa anomaly map** — Eastern_395.h5 (SGL survey, upward-continued to 395 m):
- Coverage: lat [44.34°, 45.65°], lon [-76.68°, -74.34°]
- Range: −1463 to +7447 nT
- Bilinear interpolation; returns 0 outside coverage area

**9-term TL aircraft interference** — from `results/run_20260220_221033/logs/custom_tl_coef.txt`:
```
Permanent  (c1-c3): Bx/Bt, By/Bt, Bz/Bt
Induced    (c4-c6): Bx,    By,    Bz
Quadratic  (c7-c9): Bx²/Bt, By²/Bt, Bz²/Bt
```
> **Note on coefficient magnitude:** The TL coefficients were trained on an SGL
> fixed-wing survey aircraft with large permanent magnetization (c1≈−70 kNm,
> c3≈+89 kNm). For a small x500 drone the expected interference would be
> ~100–1000 nT. Use `--no-tl` to isolate anomaly-only simulation, or retrain
> TL coefficients on a drone-specific calibration flight.

**Noise:** Gaussian Normal(0, σ) nT where σ=5.0 nT by default (matches IIS2MDC spec).

## Quickstart

### 1. Export Ottawa map (once)
```bash
julia gazebo/plugins/export_ottawa_map.jl
```
Produces `gazebo/plugins/data/ottawa_map.npz` (~240 MB).

### 2. Install Python dependencies
```bash
pip install numpy scipy pyulog
# gz-transport13 Python bindings come with Gazebo Harmonic apt package
```

### 3. Validate (no Gazebo needed)
```bash
python3 gazebo/plugins/magnetic_map_bridge.py --dry-run
```

### 4. Run full simulation
```bash
# Terminal A — bridge
python3 gazebo/plugins/magnetic_map_bridge.py --verbose

# Terminal B — PX4 SITL
export GZ_SIM_RESOURCE_PATH=$GZ_SIM_RESOURCE_PATH:$(pwd)/gazebo/models
PX4_HOME_LAT=44.76 PX4_HOME_LON=-76.15 \
  make -C ~/PX4-Autopilot px4_sitl gz_x500_magnav

# Or use the combined launcher:
bash gazebo/plugins/run_magnav_sim.sh
```

### 5. Fly box mission
Load `gazebo/scripts/box_mission.plan` in QGroundControl.

### 6. Analyse the log
```bash
# Find latest ULG:
ULG=$(ls -t ~/PX4-Autopilot/build/px4_sitl_default/logs/*.ulg | head -1)

# Run TL validation:
python3 gazebo/scripts/test_tl_on_sim.py "$ULG"
```

## Expected Results

| Metric | Before TL | After TL |
|--------|-----------|----------|
| Heading correlation | high (>0.1) | low (<0.05) |
| Variance reduction | — | >30% |
| Est. nav DRMS | ~3× GPS | ~1× GPS |

## Bridge Options

```
--map PATH       NPZ map file           (default: gazebo/plugins/data/ottawa_map.npz)
--coef PATH      TL coefficient file    (default: results/.../custom_tl_coef.txt)
--noise FLOAT    Noise std-dev (nT)     (default: 5.0)
--world STR      Gazebo world name      (default: default)
--model STR      Vehicle model name     (default: x500_magnav)
--no-tl          Disable TL interference
--no-anomaly     Disable anomaly map
--verbose        Print per-sample debug
--dry-run        Validate config, then exit
```

## Topic Reference

| Topic | Type | Direction | Description |
|-------|------|-----------|-------------|
| `.../navsat_sensor/navsat` | gz NavSat | subscribe | GPS position from Gazebo |
| `.../magnetometer_sensor/magnetometer` | gz Magnetometer | subscribe | Raw Earth field (WMM) |
| `/magnav/magnetometer` | gz Magnetometer | publish | Earth + anomaly + TL |
| `/magnav/mag_anomaly_nT` | gz Float | publish | Anomaly value at position |
| `/magnav/tl_interference_nT` | gz Float | publish | TL scalar interference |
