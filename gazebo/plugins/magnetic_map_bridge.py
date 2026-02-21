#!/usr/bin/env python3
"""
magnetic_map_bridge.py — Gazebo MagNav Magnetic Map Bridge

Stand-alone Python process that injects realistic magnetic anomaly +
aircraft TL interference into the Gazebo Harmonic magnetometer stream.

Architecture
────────────
  Gazebo SITL publishes:
    /world/.../navsat     → GPS position  (gz NavSat msg)
    /world/.../magnetometer → raw Earth field (gz Magnetometer msg)

  This bridge:
    1. Loads Ottawa anomaly map  (NPZ grid, built by export_ottawa_map.jl)
    2. Loads 9-term TL coefficients (custom_tl_coef.txt)
    3. For every mag update, computes:
         total = B_earth(raw) + B_anomaly(lat,lon) + B_tl(aircraft) + noise
    4. Publishes:
         /magnav/magnetometer    → gz Magnetometer  (field_tesla modified)
         /magnav/mag_anomaly_nT  → gz Float          (anomaly scalar only)
         /magnav/tl_interference_nT → gz Float       (TL term scalar only)

PX4 / downstream consumers subscribe to /magnav/magnetometer instead
of the default mag topic.  The test_tl_on_sim.py script then reads the
.ulg log that PX4 wrote from this modified source.

Usage
─────
  # 1. Export map (once):
  julia gazebo/plugins/export_ottawa_map.jl

  # 2. Run bridge alongside PX4 SITL:
  python3 gazebo/plugins/magnetic_map_bridge.py [options]

  # 3. Record / fly:
  PX4_HOME_LAT=44.76 PX4_HOME_LON=-76.15 make px4_sitl gz_x500

Options
───────
  --map      PATH  NPZ map file  (default: gazebo/plugins/data/ottawa_map.npz)
  --coef     PATH  TL coef file  (default: results/run_20260220_221033/logs/custom_tl_coef.txt)
  --noise    FLOAT Gaussian noise std-dev in nT  (default: 5.0 nT)
  --world    STR   Gazebo world name  (default: default)
  --model    STR   Vehicle model name (default: x500_magnav)
  --no-tl         Disable TL interference (anomaly map only)
  --no-anomaly    Disable anomaly map     (TL only)
  --verbose        Print per-sample debug
  --dry-run        Load map & coef, validate, then exit

Dependencies
────────────
  pip install numpy scipy
  gz-transport13 Python bindings (from Gazebo Harmonic apt package)
"""

import argparse
import math
import sys
import time
import threading
from pathlib import Path
from typing import Optional

import numpy as np
from scipy.interpolate import RegularGridInterpolator

# ── Gazebo Python bindings ────────────────────────────────────────────────────
try:
    from gz.transport13 import Node
    from gz.msgs10.magnetometer_pb2 import Magnetometer
    from gz.msgs10.navsat_pb2 import NavSat
    from gz.msgs10.float_pb2 import Float as GzFloat
except ImportError as e:
    sys.exit(
        f"Gazebo Python bindings not found: {e}\n"
        "Install Gazebo Harmonic and ensure gz-transport13 Python bindings are available.\n"
        "  sudo apt install gz-harmonic"
    )

# ── Repository root (two levels up from this file) ───────────────────────────
_REPO_ROOT = Path(__file__).resolve().parents[2]

DEFAULT_MAP_PATH  = _REPO_ROOT / "gazebo/plugins/data/ottawa_map.npz"
DEFAULT_COEF_PATH = _REPO_ROOT / "results/run_20260220_221033/logs/custom_tl_coef.txt"

TESLA_TO_NT = 1e9   # 1 T = 1e9 nT
NT_TO_TESLA = 1e-9

# ─────────────────────────────────────────────────────────────────────────────
# Coefficient loading
# ─────────────────────────────────────────────────────────────────────────────

def load_tl_coef(path: Path) -> np.ndarray:
    """
    Parse custom_tl_coef.txt → 9-element float64 array.

    Format (annotated or bare):
        c[01] = -69811.81361053
    """
    coef = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line.startswith("c["):
                continue
            value_part = line.split("=", 1)[1]
            coef.append(float(value_part.split()[0]))
    if len(coef) != 9:
        raise ValueError(
            f"Expected 9 TL coefficients in {path}, found {len(coef)}."
        )
    return np.array(coef, dtype=np.float64)


# ─────────────────────────────────────────────────────────────────────────────
# 9-term TL design matrix
# ─────────────────────────────────────────────────────────────────────────────

def tl_interference(Bx_nT: float, By_nT: float, Bz_nT: float,
                    coef: np.ndarray) -> float:
    """
    Compute the scalar aircraft interference (nT) from the 9-term TL model.

    Design matrix matches custom_tl.jl exactly:
      cols 0-2  permanent:         Bx/Bt, By/Bt, Bz/Bt
      cols 3-5  linear induced:    Bx, By, Bz
      cols 6-8  quadratic induced: Bx²/Bt, By²/Bt, Bz²/Bt

    Parameters
    ----------
    Bx_nT, By_nT, Bz_nT : float
        Earth background field components in nT (from Gazebo mag sensor).
    coef : np.ndarray shape (9,)
        Fitted TL coefficients from custom_tl_coef.txt.

    Returns
    -------
    float
        Predicted aircraft-induced interference in nT.
    """
    Bt = math.sqrt(Bx_nT**2 + By_nT**2 + Bz_nT**2)
    Bt_safe = max(Bt, 1.0)   # 1 nT floor (matches Julia code)

    a = np.array([
        Bx_nT / Bt_safe,          # permanent
        By_nT / Bt_safe,
        Bz_nT / Bt_safe,
        Bx_nT,                     # linear induced
        By_nT,
        Bz_nT,
        Bx_nT**2 / Bt_safe,        # quadratic induced
        By_nT**2 / Bt_safe,
        Bz_nT**2 / Bt_safe,
    ], dtype=np.float64)

    return float(a @ coef)


# ─────────────────────────────────────────────────────────────────────────────
# Ottawa anomaly map
# ─────────────────────────────────────────────────────────────────────────────

class OttawaMap:
    """
    Loads the Ottawa magnetic anomaly map from NPZ and provides
    fast bilinear interpolation at arbitrary (lat, lon) positions.

    The NPZ file is produced by export_ottawa_map.jl and contains:
        lat_rad, lon_rad — 1D grid vectors (radians)
        map_nT           — 2D anomaly array  (lat × lon), NaN outside survey
        alt_m            — survey altitude (scalar)

    Query method
    ────────────
    The gz NavSat message gives latitude/longitude in *degrees*.  This class
    accepts degrees and converts internally.
    """

    def __init__(self, npz_path):
        npz_path = Path(npz_path)
        if not npz_path.exists():
            raise FileNotFoundError(
                f"Ottawa map NPZ not found: {npz_path}\n"
                "Run first:  julia gazebo/plugins/export_ottawa_map.jl"
            )
        data = np.load(str(npz_path))
        self.lat_rad  = data["lat_rad"]          # (N_lat,)
        self.lon_rad  = data["lon_rad"]          # (N_lon,)
        self.map_nT   = data["map_nT"]           # (N_lat, N_lon)
        self.alt_m    = float(data["alt_m"][0])

        # Build interpolator — RegularGridInterpolator is O(log n) per query
        # bounds_error=False + fill_value=0 → returns 0 outside coverage area
        self._itp = RegularGridInterpolator(
            (self.lat_rad, self.lon_rad),
            self.map_nT,
            method="linear",
            bounds_error=False,
            fill_value=0.0,         # 0 anomaly outside map coverage
        )

        # Coverage in degrees for display
        self.lat_min_deg = math.degrees(self.lat_rad[0])
        self.lat_max_deg = math.degrees(self.lat_rad[-1])
        self.lon_min_deg = math.degrees(self.lon_rad[0])
        self.lon_max_deg = math.degrees(self.lon_rad[-1])

    def query(self, lat_deg: float, lon_deg: float) -> float:
        """
        Return anomaly value at (lat_deg, lon_deg) in nT.
        Returns 0.0 outside survey coverage.
        """
        lat_r = math.radians(lat_deg)
        lon_r = math.radians(lon_deg)
        # RegularGridInterpolator expects a (1, 2) array for a single point
        val = self._itp([[lat_r, lon_r]])[0]
        # Replace NaN (masked cells) with 0
        return 0.0 if math.isnan(val) else float(val)

    def in_coverage(self, lat_deg: float, lon_deg: float) -> bool:
        lat_r = math.radians(lat_deg)
        lon_r = math.radians(lon_deg)
        return (self.lat_rad[0] <= lat_r <= self.lat_rad[-1] and
                self.lon_rad[0] <= lon_r <= self.lon_rad[-1])

    def __repr__(self):
        return (
            f"OttawaMap("
            f"lat=[{self.lat_min_deg:.3f},{self.lat_max_deg:.3f}]°, "
            f"lon=[{self.lon_min_deg:.3f},{self.lon_max_deg:.3f}]°, "
            f"alt={self.alt_m:.0f}m, "
            f"shape={self.map_nT.shape})"
        )


# ─────────────────────────────────────────────────────────────────────────────
# Bridge state (updated by subscriber callbacks, read by publisher timer)
# ─────────────────────────────────────────────────────────────────────────────

class BridgeState:
    """Thread-safe shared state between subscriber callbacks and publisher."""

    def __init__(self):
        self._lock = threading.Lock()
        # GPS
        self.lat_deg: float = 44.76   # Ottawa default
        self.lon_deg: float = -76.15
        self.alt_m:   float = 395.0
        self.gps_valid: bool = False

        # Magnetometer (Tesla, from gz sensor)
        self.Bx_T: float = 0.0
        self.By_T: float = 0.0
        self.Bz_T: float = 0.0
        self.mag_valid: bool = False

        # Counters
        self.n_gps: int = 0
        self.n_mag: int = 0
        self.n_pub: int = 0

    def update_gps(self, lat: float, lon: float, alt: float):
        with self._lock:
            self.lat_deg = lat
            self.lon_deg = lon
            self.alt_m   = alt
            self.gps_valid = True
            self.n_gps += 1

    def update_mag(self, bx: float, by: float, bz: float):
        with self._lock:
            self.Bx_T = bx
            self.By_T = by
            self.Bz_T = bz
            self.mag_valid = True
            self.n_mag += 1

    def snapshot(self):
        with self._lock:
            return (
                self.lat_deg, self.lon_deg, self.alt_m,
                self.Bx_T, self.By_T, self.Bz_T,
                self.gps_valid, self.mag_valid,
            )

    def inc_pub(self):
        with self._lock:
            self.n_pub += 1

    def stats(self):
        with self._lock:
            return self.n_gps, self.n_mag, self.n_pub


# ─────────────────────────────────────────────────────────────────────────────
# Topic helpers
# ─────────────────────────────────────────────────────────────────────────────

def make_sensor_topic(world: str, model: str, sensor_type: str,
                      sensor_name: str) -> str:
    """
    Construct the canonical gz-transport topic for a Gazebo Harmonic sensor.

    Pattern:
      /world/{world}/model/{model}/link/base_link/sensor/{sensor_name}/{sensor_type}

    Gazebo Harmonic (gz-sim8) default topic naming convention.
    """
    return (
        f"/world/{world}/model/{model}/link/base_link"
        f"/sensor/{sensor_name}/{sensor_type}"
    )


# ─────────────────────────────────────────────────────────────────────────────
# Main bridge class
# ─────────────────────────────────────────────────────────────────────────────

class MagneticMapBridge:
    """
    Subscribes to Gazebo GPS + magnetometer topics, injects Ottawa anomaly
    and TL interference, and republishes on /magnav/* topics.

    Signal model
    ────────────
    The gz magnetometer sensor gives the Earth's background field as a 3D
    vector [Bx, By, Bz] in Tesla (WMM model).

    We compute a modified *scalar* total-field measurement:

        Bt_earth   = sqrt(Bx² + By² + Bz²)          [nT]
        Bt_anomaly = OttawaMap.query(lat, lon)        [nT]
        Bt_tl      = tl_interference(Bx, By, Bz, c)  [nT]
        noise      = Normal(0, sigma)                 [nT]

        Bt_total = Bt_earth + Bt_anomaly + Bt_tl + noise

    The published /magnav/magnetometer Magnetometer message contains a
    field_tesla vector scaled so that its magnitude equals Bt_total, with
    the same direction as the original Earth field vector.  This is the
    convention used by the PX4 mag sensor (magnitude = total field strength).

    Note: A real fluxgate magnetometer measures the vector field.  The TL
    interference is a scalar addition to the total intensity, which is the
    standard MagNav formulation (as in MagNav.jl / custom_tl.jl).
    """

    def __init__(
        self,
        ottawa_map: OttawaMap,
        tl_coef: Optional[np.ndarray],
        noise_std_nT: float = 5.0,
        world: str = "default",
        model: str = "x500_magnav",
        use_tl: bool = True,
        use_anomaly: bool = True,
        verbose: bool = False,
    ):
        self.map         = ottawa_map
        self.coef        = tl_coef
        self.noise_std   = noise_std_nT
        self.use_tl      = use_tl and (tl_coef is not None)
        self.use_anomaly = use_anomaly
        self.verbose     = verbose
        self.rng         = np.random.default_rng(seed=42)

        self.state = BridgeState()

        # ── gz-transport node ─────────────────────────────────────────────────
        self.node = Node()

        # ── Subscribe: NavSat (GPS) ───────────────────────────────────────────
        navsat_topic = make_sensor_topic(world, model, "navsat", "navsat_sensor")
        ok = self.node.subscribe(NavSat, navsat_topic, self._on_navsat)
        if not ok:
            print(f"[WARN] Could not subscribe to NavSat topic: {navsat_topic}")
            print("       Will use default Ottawa position until GPS is received.")
        else:
            print(f"[INFO] Subscribed to GPS:  {navsat_topic}")

        # ── Subscribe: Magnetometer ───────────────────────────────────────────
        mag_topic = make_sensor_topic(world, model, "magnetometer",
                                      "magnetometer_sensor")
        ok = self.node.subscribe(Magnetometer, mag_topic, self._on_magnetometer)
        if not ok:
            print(f"[WARN] Could not subscribe to Magnetometer topic: {mag_topic}")
        else:
            print(f"[INFO] Subscribed to Mag:  {mag_topic}")

        # ── Advertise output topics ───────────────────────────────────────────
        self._pub_mag    = self.node.advertise("/magnav/magnetometer",
                                                Magnetometer)
        self._pub_anomaly = self.node.advertise("/magnav/mag_anomaly_nT", GzFloat)
        self._pub_tl      = self.node.advertise("/magnav/tl_interference_nT",
                                                GzFloat)

        print(f"[INFO] Publishing to: /magnav/magnetometer")
        print(f"[INFO] Publishing to: /magnav/mag_anomaly_nT")
        print(f"[INFO] Publishing to: /magnav/tl_interference_nT")

    # ── Subscriber callbacks ──────────────────────────────────────────────────

    def _on_navsat(self, msg: NavSat):
        self.state.update_gps(msg.latitude_deg, msg.longitude_deg, msg.altitude)

    def _on_magnetometer(self, msg: Magnetometer):
        self.state.update_mag(
            msg.field_tesla.x,
            msg.field_tesla.y,
            msg.field_tesla.z,
        )
        # Each incoming mag message triggers one output publication
        self._publish_once()

    # ── Core computation ──────────────────────────────────────────────────────

    def _publish_once(self):
        lat, lon, alt, Bx_T, By_T, Bz_T, gps_ok, mag_ok = self.state.snapshot()

        if not mag_ok:
            return   # nothing to publish yet

        # Convert to nT for all calculations
        Bx_nT = Bx_T * TESLA_TO_NT
        By_nT = By_T * TESLA_TO_NT
        Bz_nT = Bz_T * TESLA_TO_NT
        Bt_earth_nT = math.sqrt(Bx_nT**2 + By_nT**2 + Bz_nT**2)

        # 1) Ottawa anomaly map value
        anomaly_nT = 0.0
        if self.use_anomaly:
            anomaly_nT = self.map.query(lat, lon)

        # 2) TL aircraft interference (9-term model)
        tl_nT = 0.0
        if self.use_tl and self.coef is not None:
            tl_nT = tl_interference(Bx_nT, By_nT, Bz_nT, self.coef)

        # 3) Gaussian sensor noise
        noise_nT = float(self.rng.normal(0.0, self.noise_std))

        # 4) Total scalar field
        Bt_total_nT = Bt_earth_nT + anomaly_nT + tl_nT + noise_nT

        # Scale the Earth-field direction vector so its magnitude = Bt_total
        # (preserves direction of B-field, adjusts magnitude only)
        if Bt_earth_nT > 1.0:   # avoid divide-by-zero
            scale = Bt_total_nT / Bt_earth_nT
        else:
            scale = 1.0

        Bx_out_T = Bx_T * scale
        By_out_T = By_T * scale
        Bz_out_T = Bz_T * scale

        # ── Publish /magnav/magnetometer ─────────────────────────────────────
        out_mag = Magnetometer()
        out_mag.field_tesla.x = Bx_out_T
        out_mag.field_tesla.y = By_out_T
        out_mag.field_tesla.z = Bz_out_T
        self._pub_mag.publish(out_mag)

        # ── Publish anomaly scalar ────────────────────────────────────────────
        out_a = GzFloat()
        out_a.data = anomaly_nT
        self._pub_anomaly.publish(out_a)

        # ── Publish TL scalar ─────────────────────────────────────────────────
        out_tl = GzFloat()
        out_tl.data = tl_nT
        self._pub_tl.publish(out_tl)

        self.state.inc_pub()

        if self.verbose:
            in_cov = self.map.in_coverage(lat, lon)
            print(
                f"[MAG] lat={lat:.5f}° lon={lon:.5f}° "
                f"{'IN' if in_cov else 'OUT'} | "
                f"earth={Bt_earth_nT:.1f} "
                f"anom={anomaly_nT:+.1f} "
                f"tl={tl_nT:+.1f} "
                f"noise={noise_nT:+.1f} "
                f"total={Bt_total_nT:.1f} nT"
            )

    # ── Status printer ────────────────────────────────────────────────────────

    def print_status(self, interval_s: float = 10.0):
        """Background thread: print stats every interval_s seconds."""
        while True:
            time.sleep(interval_s)
            n_gps, n_mag, n_pub = self.state.stats()
            lat, lon, alt, *_ = self.state.snapshot()
            in_cov = self.map.in_coverage(lat, lon)
            anomaly = self.map.query(lat, lon)
            print(
                f"[STATUS] GPS msgs={n_gps} | Mag msgs={n_mag} | "
                f"Published={n_pub} | "
                f"pos=({lat:.4f}°,{lon:.4f}°) "
                f"{'[IN MAP]' if in_cov else '[OUT OF MAP - anomaly=0]'} | "
                f"anomaly={anomaly:+.1f} nT"
            )


# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Gazebo MagNav magnetic map bridge",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--map",   default=str(DEFAULT_MAP_PATH),
                   help=f"NPZ map file (default: {DEFAULT_MAP_PATH})")
    p.add_argument("--coef",  default=str(DEFAULT_COEF_PATH),
                   help=f"TL coefficient file (default: {DEFAULT_COEF_PATH})")
    p.add_argument("--noise", type=float, default=5.0,
                   help="Gaussian noise std-dev in nT (default: 5.0)")
    p.add_argument("--world", default="default",
                   help="Gazebo world name (default: default)")
    p.add_argument("--model", default="x500_magnav",
                   help="Vehicle model name (default: x500_magnav)")
    p.add_argument("--no-tl",      action="store_true",
                   help="Disable TL interference (anomaly map only)")
    p.add_argument("--no-anomaly", action="store_true",
                   help="Disable anomaly map (TL only)")
    p.add_argument("--verbose",    action="store_true",
                   help="Print per-sample debug output")
    p.add_argument("--dry-run",    action="store_true",
                   help="Load map & coef, validate, then exit")
    return p.parse_args()


def main():
    args = parse_args()

    print("=" * 60)
    print("  Gazebo MagNav Magnetic Map Bridge")
    print("=" * 60)

    # ── Load Ottawa map ───────────────────────────────────────────────────────
    map_path = Path(args.map)
    print(f"\n▶ Loading Ottawa anomaly map: {map_path}")
    try:
        ottawa = OttawaMap(map_path)
    except FileNotFoundError as e:
        sys.exit(str(e))
    print(f"  {ottawa}")

    # Quick sanity check: Ottawa city centre
    test_lat, test_lon = 44.76, -76.15
    test_val = ottawa.query(test_lat, test_lon)
    in_cov   = ottawa.in_coverage(test_lat, test_lon)
    print(f"  Spot check Ottawa ({test_lat}°, {test_lon}°): "
          f"{test_val:.1f} nT  [{'in' if in_cov else 'OUT OF'} coverage]")

    # ── Load TL coefficients ──────────────────────────────────────────────────
    coef_path = Path(args.coef)
    coef: Optional[np.ndarray] = None
    if args.no_tl:
        print("\n▶ TL interference: DISABLED (--no-tl)")
    else:
        print(f"\n▶ Loading TL coefficients: {coef_path}")
        try:
            coef = load_tl_coef(coef_path)
            print(f"  Coefficients (9-term): {np.round(coef, 4)}")
        except (FileNotFoundError, ValueError) as e:
            sys.exit(f"TL coefficient load failed: {e}")

    if args.no_anomaly:
        print("▶ Anomaly map: DISABLED (--no-anomaly)")

    print(f"\n▶ Noise std-dev: {args.noise:.1f} nT")

    if args.dry_run:
        print("\n✓ Dry-run complete — map and coefficients loaded OK.")
        return

    # ── Create bridge ─────────────────────────────────────────────────────────
    print(f"\n▶ Starting bridge (world='{args.world}', model='{args.model}')")
    bridge = MagneticMapBridge(
        ottawa_map   = ottawa,
        tl_coef      = coef,
        noise_std_nT = args.noise,
        world        = args.world,
        model        = args.model,
        use_tl       = not args.no_tl,
        use_anomaly  = not args.no_anomaly,
        verbose      = args.verbose,
    )

    # ── Background status thread ──────────────────────────────────────────────
    status_thread = threading.Thread(
        target=bridge.print_status, args=(10.0,), daemon=True
    )
    status_thread.start()

    print("\n✓ Bridge running.  Press Ctrl-C to stop.\n")
    print("  Waiting for Gazebo to publish sensor data...")
    print("  (Start PX4 SITL if not already running)\n")

    try:
        while True:
            time.sleep(0.5)
    except KeyboardInterrupt:
        n_gps, n_mag, n_pub = bridge.state.stats()
        print(f"\n[INFO] Shutting down.  "
              f"GPS={n_gps}  Mag={n_mag}  Published={n_pub}")


if __name__ == "__main__":
    main()
