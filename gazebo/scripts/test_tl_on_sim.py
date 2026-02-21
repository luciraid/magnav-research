#!/usr/bin/env python3
"""
test_tl_on_sim.py — Test the trained 9-term TL model on a PX4 Gazebo flight log.

Applies saved Tolles-Lawson coefficients to a .ulg log, computes magnetic
compensation quality, and generates a 4-panel plot comparing GPS-only and
TL-compensated navigation estimates.

Usage:
    python3 gazebo/scripts/test_tl_on_sim.py <path/to/flight.ulg>
    python3 gazebo/scripts/test_tl_on_sim.py flight.ulg --output gazebo/results/
    python3 gazebo/scripts/test_tl_on_sim.py flight.ulg \\
        --coef results/run_20260220_221033/logs/custom_tl_coef.txt

What is computed (pure Python, no Julia, no MagNav):
    - Magnetic scalar before/after TL compensation
    - Heading-correlated interference (root cause of navigation error in MagNav)
    - Estimated DRMS improvement, extrapolated from the real Flt1006 benchmark
    - GPS trajectory coloured by altitude

What requires MagNav.jl + anomaly map for exact results:
    - True EKF navigation DRMS (run benchmark/bench_tl_compare.jl on real data)

Dependencies:
    pip install pyulog numpy matplotlib scipy
"""

import sys
import os
import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")           # headless — no display needed
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
from scipy.interpolate import interp1d

try:
    from pyulog import ULog
except ImportError:
    sys.exit("pyulog not found.  Install with:  pip install pyulog")


# ── Constants ─────────────────────────────────────────────────────────────────

GAUSS_TO_NT   = 1e5          # 1 Gauss = 100 000 nT
TARGET_HZ     = 10           # resample rate (matches SGL dataset)
GPS_DRMS_M    = 3.0          # Gazebo SITL GPS noise baseline (m)

# Reference values from real Flt1006 benchmark (run_20260220_221033)
REF_DRMS_OFFICIAL  = 0.59    # m  (official MagNav TL)
REF_DRMS_CUSTOM    = 0.36    # m  (custom 9-term TL)
REF_HDGCORR_OFFIC  = 0.027   # |corr| after official TL (from README)
REF_HDGCORR_CUSTOM = 0.015   # |corr| after custom TL   (from README)

DEFAULT_COEF = Path(__file__).resolve().parents[2] / \
               "results/run_20260220_221033/logs/custom_tl_coef.txt"


# ── Coefficient parser ────────────────────────────────────────────────────────

def load_coef(path: Path) -> np.ndarray:
    """
    Parse custom_tl_coef.txt and return a 9-element float64 array.

    Handles both the annotated format (with trailing comment):
        c[01] = -69811.81361053    (Permanent: ...)
    and the bare format:
        c[01] = -69811.81361053
    """
    coef = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line.startswith("c["):
                continue
            # Take only the part after '='
            value_part = line.split("=", 1)[1]
            # Strip trailing comment if present
            numeric = value_part.split()[0]
            coef.append(float(numeric))
    if len(coef) != 9:
        raise ValueError(
            f"Expected 9 TL coefficients in {path}, found {len(coef)}.\n"
            "Check that the file matches the custom_tl_coef.txt format."
        )
    return np.array(coef, dtype=np.float64)


# ── 9-term TL design matrix ───────────────────────────────────────────────────

def build_A9(Bx: np.ndarray, By: np.ndarray, Bz: np.ndarray) -> np.ndarray:
    """
    Build the 9-term Tolles-Lawson design matrix (n × 9).

    Matches custom_tl.jl exactly:
      cols 0-2  permanent:         Bx/Bt, By/Bt, Bz/Bt
      cols 3-5  linear induced:    Bx, By, Bz
      cols 6-8  quadratic induced: Bx²/Bt, By²/Bt, Bz²/Bt

    All in nT (inputs must be in nT).
    """
    Bt      = np.sqrt(Bx**2 + By**2 + Bz**2)
    Bt_safe = np.maximum(Bt, 1.0)       # 1 nT floor, matches Julia code

    A = np.column_stack([
        Bx / Bt_safe,           # permanent
        By / Bt_safe,
        Bz / Bt_safe,
        Bx,                     # linear induced
        By,
        Bz,
        Bx**2 / Bt_safe,        # quadratic induced
        By**2 / Bt_safe,
        Bz**2 / Bt_safe,
    ])
    return A


# ── Signal helpers ────────────────────────────────────────────────────────────

def resample(t_src: np.ndarray, y_src: np.ndarray, t_dst: np.ndarray) -> np.ndarray:
    """Linear interpolation of y_src onto t_dst."""
    f = interp1d(t_src, y_src, kind="linear",
                 bounds_error=False, fill_value=(y_src[0], y_src[-1]))
    return f(t_dst)


def detrend_linear(y: np.ndarray) -> np.ndarray:
    """Remove a least-squares linear trend from y."""
    n = len(y)
    t = np.arange(n, dtype=float)
    A = np.column_stack([np.ones(n), t])
    c, _, _, _ = np.linalg.lstsq(A, y, rcond=None)
    return y - A @ c


def heading_correlation(mag: np.ndarray, heading_rad: np.ndarray) -> float:
    """|Pearson correlation| between mag signal and heading."""
    return float(np.abs(np.corrcoef(mag, heading_rad)[0, 1]))


def rms(x: np.ndarray) -> float:
    return float(np.sqrt(np.mean(x**2)))


# ── ULog extraction ───────────────────────────────────────────────────────────

def quat_to_yaw(q: np.ndarray) -> np.ndarray:
    """
    Convert PX4 quaternion array (n × 4, columns w x y z) to yaw (radians).
    Yaw = atan2(2(wz + xy), 1 - 2(y² + z²))
    """
    w, x, y, z = q[:, 0], q[:, 1], q[:, 2], q[:, 3]
    return np.arctan2(2*(w*z + x*y), 1 - 2*(y**2 + z**2))


def extract_ulg(ulg_path: Path) -> dict:
    """
    Extract GPS, fluxgate magnetometer, and heading from a PX4 .ulg file.

    Returns a dict with keys:
        t_s, lat_deg, lon_deg, alt_m,
        Bx_nT, By_nT, Bz_nT,
        heading_rad
    All arrays resampled to TARGET_HZ on a common time grid.

    Raises RuntimeError if required datasets are missing.
    """
    log = ULog(str(ulg_path))
    available = {d.name for d in log.data_list}

    # ── GPS ───────────────────────────────────────────────────────────────────
    if "vehicle_gps_position" not in available:
        raise RuntimeError("Missing 'vehicle_gps_position' dataset in ULog.")
    gps = log.get_dataset("vehicle_gps_position").data
    t_gps  = gps["timestamp"] / 1e6   # μs → s
    lat    = gps["latitude_deg"]       # already degrees
    lon    = gps["longitude_deg"]
    alt    = gps["altitude_msl_m"]     # already metres

    # ── Magnetometer ──────────────────────────────────────────────────────────
    mag_topic = "vehicle_magnetometer" if "vehicle_magnetometer" in available \
                else "sensor_mag"
    if mag_topic not in available:
        raise RuntimeError(
            "No magnetometer dataset found. "
            "Expected 'vehicle_magnetometer' or 'sensor_mag'."
        )
    mag    = log.get_dataset(mag_topic).data
    t_mag  = mag["timestamp"] / 1e6

    # Field names vary between PX4 versions
    if "magnetometer_ga[0]" in mag:
        mx = mag["magnetometer_ga[0]"] * GAUSS_TO_NT
        my = mag["magnetometer_ga[1]"] * GAUSS_TO_NT
        mz = mag["magnetometer_ga[2]"] * GAUSS_TO_NT
    elif "x" in mag:
        # Older PX4: values in Gauss
        mx = mag["x"] * GAUSS_TO_NT
        my = mag["y"] * GAUSS_TO_NT
        mz = mag["z"] * GAUSS_TO_NT
    else:
        raise RuntimeError(f"Cannot find magnetometer field names in '{mag_topic}'.")

    # ── Attitude / heading ────────────────────────────────────────────────────
    if "vehicle_attitude" not in available:
        raise RuntimeError("Missing 'vehicle_attitude' dataset in ULog.")
    att   = log.get_dataset("vehicle_attitude").data
    t_att = att["timestamp"] / 1e6

    if "q[0]" in att:
        # PX4 stores quaternion (w, x, y, z)
        q = np.column_stack([att["q[0]"], att["q[1]"], att["q[2]"], att["q[3]"]])
        yaw_raw = quat_to_yaw(q)
    elif "yaw" in att:
        yaw_raw = att["yaw"]
    else:
        raise RuntimeError("Cannot extract heading from 'vehicle_attitude'.")

    # ── Build common time grid ────────────────────────────────────────────────
    t0     = max(t_gps[0], t_mag[0], t_att[0])
    t1     = min(t_gps[-1], t_mag[-1], t_att[-1])
    dt     = 1.0 / TARGET_HZ
    t_grid = np.arange(t0, t1, dt)

    if len(t_grid) < 10:
        raise RuntimeError(
            f"Flight too short after sync ({len(t_grid)} samples). "
            "Check that all datasets overlap in time."
        )

    return {
        "t_s":        t_grid - t_grid[0],
        "lat_deg":    resample(t_gps, lat,      t_grid),
        "lon_deg":    resample(t_gps, lon,      t_grid),
        "alt_m":      resample(t_gps, alt,      t_grid),
        "Bx_nT":      resample(t_mag, mx,       t_grid),
        "By_nT":      resample(t_mag, my,       t_grid),
        "Bz_nT":      resample(t_mag, mz,       t_grid),
        "heading_rad":resample(t_att, yaw_raw,  t_grid),
    }


# ── DRMS estimate ─────────────────────────────────────────────────────────────

def estimate_drms(hdg_corr: float, reference_drms: float,
                  reference_hdg_corr: float) -> float:
    """
    Extrapolate an expected navigation DRMS from a heading-correlation measurement.

    Physical basis: In MagNav EKF navigation, DRMS scales approximately with
    the heading-correlated component of the magnetic error (empirically confirmed
    by the ratio 0.36/0.59 ≈ 0.015/0.027 from the Flt1006 benchmark).

    This estimate assumes:
      - Same magnetic anomaly map quality as Flt1006
      - Same EKF parameters (R=100, same P0/Qd)
    Label as ESTIMATED on all plots.

    Returns estimated DRMS in metres.
    """
    if reference_hdg_corr == 0:
        return 0.0
    return reference_drms * (hdg_corr / reference_hdg_corr)


# ── Plotting ──────────────────────────────────────────────────────────────────

def make_plot(data: dict, Bt_raw: np.ndarray, Bt_comp: np.ndarray,
              interference: np.ndarray, metrics: dict,
              ulg_name: str, output_dir: Path) -> Path:
    """
    Generate a 4-panel summary figure.

    Panel 1 (top-left):  Magnetic scalar — raw vs TL-compensated over time
    Panel 2 (top-right): GPS trajectory coloured by altitude
    Panel 3 (bot-left):  Heading vs detrended B scatter (shows interference pattern)
    Panel 4 (bot-right): DRMS bar chart (GPS / raw-mag / TL-compensated estimates)
    """
    t      = data["t_s"]
    hdg    = np.degrees(data["heading_rad"]) % 360
    lat    = data["lat_deg"]
    lon    = data["lon_deg"]
    alt    = data["alt_m"]

    Bt_dt_raw  = detrend_linear(Bt_raw)
    Bt_dt_comp = detrend_linear(Bt_comp)

    fig = plt.figure(figsize=(14, 10))
    gs  = gridspec.GridSpec(2, 2, figure=fig,
                             hspace=0.38, wspace=0.30,
                             left=0.08, right=0.97,
                             top=0.92, bottom=0.08)

    # ── Panel 1: Magnetic scalar over time ───────────────────────────────────
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.plot(t / 60, Bt_raw,  color="#3a86ff", lw=0.8, alpha=0.7, label="Raw |B|")
    ax1.plot(t / 60, Bt_comp, color="#ff6b35", lw=1.2, label="TL-compensated |B|")
    ax1.set_xlabel("Time (min)")
    ax1.set_ylabel("Scalar field (nT)")
    ax1.set_title("Magnetic Scalar — Raw vs TL-Compensated")
    ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.3)

    # Add variance reduction annotation
    var_red = metrics["var_reduction_pct"]
    ax1.text(0.02, 0.96, f"Variance reduction: {var_red:.1f}%",
             transform=ax1.transAxes, fontsize=8,
             verticalalignment="top",
             bbox=dict(boxstyle="round,pad=0.3", facecolor="wheat", alpha=0.7))

    # ── Panel 2: GPS trajectory ───────────────────────────────────────────────
    ax2 = fig.add_subplot(gs[0, 1])
    sc = ax2.scatter(lon, lat, c=alt, cmap="viridis", s=2, zorder=3)
    plt.colorbar(sc, ax=ax2, label="Altitude (m)", pad=0.01)
    ax2.plot(lon[0],  lat[0],  "g^", ms=8, zorder=5, label="Start")
    ax2.plot(lon[-1], lat[-1], "rs", ms=8, zorder=5, label="End")
    ax2.set_xlabel("Longitude (°)")
    ax2.set_ylabel("Latitude (°)")
    ax2.set_title("GPS Trajectory")
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)
    duration_min = t[-1] / 60
    ax2.text(0.02, 0.04,
             f"{len(t)} samples  |  {duration_min:.1f} min",
             transform=ax2.transAxes, fontsize=8,
             verticalalignment="bottom",
             bbox=dict(boxstyle="round,pad=0.3", facecolor="lightcyan", alpha=0.7))

    # ── Panel 3: Heading-dependent interference ───────────────────────────────
    ax3 = fig.add_subplot(gs[1, 0])
    # Subsample if flight is long (avoid overplotting)
    step = max(1, len(t) // 2000)
    ax3.scatter(hdg[::step], Bt_dt_raw[::step],
                s=2, alpha=0.4, color="#3a86ff", label="Raw (detrended)")
    ax3.scatter(hdg[::step], Bt_dt_comp[::step],
                s=2, alpha=0.4, color="#ff6b35", label="TL-compensated")
    ax3.axhline(0, color="gray", lw=0.8, ls="--")
    ax3.set_xlabel("Aircraft heading (°)")
    ax3.set_ylabel("Detrended B (nT)")
    ax3.set_title("Heading-Correlated Interference")
    ax3.set_xlim(0, 360)
    ax3.set_xticks([0, 90, 180, 270, 360])
    ax3.legend(fontsize=8)
    ax3.grid(True, alpha=0.3)

    hc_raw  = metrics["hdg_corr_raw"]
    hc_comp = metrics["hdg_corr_comp"]
    ax3.text(0.02, 0.96,
             f"|corr| before: {hc_raw:.4f}   after: {hc_comp:.4f}",
             transform=ax3.transAxes, fontsize=8,
             verticalalignment="top",
             bbox=dict(boxstyle="round,pad=0.3", facecolor="wheat", alpha=0.7))

    # ── Panel 4: Estimated DRMS bar chart ────────────────────────────────────
    ax4 = fig.add_subplot(gs[1, 1])

    drms_gps   = GPS_DRMS_M
    drms_raw   = metrics["drms_est_raw_m"]
    drms_tl    = metrics["drms_est_tl_m"]

    bars  = ["GPS\n(baseline)", "Mag nav\n(no TL)*", "Mag nav\n(TL)*"]
    vals  = [drms_gps, drms_raw, drms_tl]
    colors = ["#06a77d", "#3a86ff", "#ff6b35"]

    rects = ax4.bar(bars, vals, color=colors, width=0.5, edgecolor="white", lw=1.2)
    for rect, val in zip(rects, vals):
        ax4.text(rect.get_x() + rect.get_width() / 2,
                 rect.get_height() + 0.01 * max(vals),
                 f"{val:.2f} m", ha="center", va="bottom", fontsize=9, fontweight="bold")

    ax4.set_ylabel("DRMS (m)")
    ax4.set_title("Navigation DRMS Comparison")
    ax4.set_ylim(0, max(vals) * 1.20)
    ax4.grid(True, axis="y", alpha=0.3)
    ax4.text(0.02, 0.01,
             "* Estimated from heading-correlation\n  scaling vs Flt1006 benchmark.\n"
             "  Exact values require MagNav.jl + EKF.",
             transform=ax4.transAxes, fontsize=7, color="dimgray",
             verticalalignment="bottom",
             bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8))

    # ── Figure title ──────────────────────────────────────────────────────────
    improvement = (drms_raw - drms_tl) / drms_raw * 100 if drms_raw > 0 else 0
    fig.suptitle(
        f"TL Compensation Validation — {ulg_name}\n"
        f"Var. reduction: {var_red:.1f}%   "
        f"Hdg corr: {hc_raw:.4f} → {hc_comp:.4f}   "
        f"Est. DRMS improvement: {improvement:.0f}%",
        fontsize=11, fontweight="bold"
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    stem     = Path(ulg_name).stem
    out_path = output_dir / f"{stem}_tl_validation.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Test trained TL model on a PX4 Gazebo .ulg flight log."
    )
    parser.add_argument("ulg", help="Path to PX4 .ulg log file")
    parser.add_argument(
        "--coef", default=str(DEFAULT_COEF),
        help=f"Path to custom_tl_coef.txt  (default: {DEFAULT_COEF})"
    )
    parser.add_argument(
        "--output", default="gazebo/results",
        help="Output directory for the plot PNG  (default: gazebo/results)"
    )
    args = parser.parse_args()

    ulg_path  = Path(args.ulg)
    coef_path = Path(args.coef)
    out_dir   = Path(args.output)

    print("=" * 60)
    print("  test_tl_on_sim.py — TL Validation on Gazebo Flight")
    print("=" * 60)

    # ── Load coefficients ─────────────────────────────────────────────────────
    print(f"\n▶ Loading TL coefficients from: {coef_path}")
    if not coef_path.exists():
        sys.exit(f"Coefficient file not found: {coef_path}")
    coef = load_coef(coef_path)
    print(f"  9 coefficients loaded: {np.round(coef, 4)}")

    # ── Extract ULog ─────────────────────────────────────────────────────────
    print(f"\n▶ Extracting data from: {ulg_path}")
    if not ulg_path.exists():
        sys.exit(f"ULog file not found: {ulg_path}")
    try:
        data = extract_ulg(ulg_path)
    except RuntimeError as e:
        sys.exit(f"ULog extraction failed: {e}")

    n = len(data["t_s"])
    dur = data["t_s"][-1]
    print(f"  Samples: {n}  ({dur/60:.1f} min at {TARGET_HZ} Hz)")
    print(f"  Lat: {data['lat_deg'].min():.4f}° – {data['lat_deg'].max():.4f}°")
    print(f"  Alt: {data['alt_m'].min():.0f} – {data['alt_m'].max():.0f} m")

    # ── Apply TL compensation ────────────────────────────────────────────────
    print("\n▶ Applying 9-term TL compensation...")
    Bx, By, Bz = data["Bx_nT"], data["By_nT"], data["Bz_nT"]
    Bt_raw      = np.sqrt(Bx**2 + By**2 + Bz**2)

    A9            = build_A9(Bx, By, Bz)
    interference  = A9 @ coef           # predicted aircraft field (nT)
    Bt_comp       = Bt_raw - interference

    print(f"  RMS interference removed: {rms(interference):.2f} nT")
    print(f"  |B| raw  — mean: {Bt_raw.mean():.1f} nT, std: {Bt_raw.std():.1f} nT")
    print(f"  |B| comp — mean: {Bt_comp.mean():.1f} nT, std: {Bt_comp.std():.1f} nT")

    # ── Metrics ───────────────────────────────────────────────────────────────
    heading = data["heading_rad"]
    Bt_dt_raw  = detrend_linear(Bt_raw)
    Bt_dt_comp = detrend_linear(Bt_comp)

    hdg_corr_raw  = heading_correlation(Bt_dt_raw,  heading)
    hdg_corr_comp = heading_correlation(Bt_dt_comp, heading)

    var_red = (1 - np.var(Bt_dt_comp) / np.var(Bt_dt_raw)) * 100

    # DRMS estimates (extrapolated from Flt1006 benchmark via heading-corr scaling)
    drms_raw_est = estimate_drms(hdg_corr_raw,  REF_DRMS_OFFICIAL,  REF_HDGCORR_OFFIC)
    drms_tl_est  = estimate_drms(hdg_corr_comp, REF_DRMS_CUSTOM,    REF_HDGCORR_CUSTOM)

    metrics = {
        "var_reduction_pct": var_red,
        "hdg_corr_raw":      hdg_corr_raw,
        "hdg_corr_comp":     hdg_corr_comp,
        "rms_interference":  rms(interference),
        "drms_est_raw_m":    drms_raw_est,
        "drms_est_tl_m":     drms_tl_est,
    }

    # ── Print summary ─────────────────────────────────────────────────────────
    print("\n── Compensation Metrics ────────────────────────────────")
    print(f"  Variance reduction        : {var_red:.1f}%")
    print(f"  Heading corr before TL    : {hdg_corr_raw:.4f}")
    print(f"  Heading corr after  TL    : {hdg_corr_comp:.4f}")
    print(f"  RMS interference removed  : {rms(interference):.2f} nT")
    print("── Estimated DRMS (extrapolated from Flt1006 benchmark) ─")
    print(f"  GPS baseline              : {GPS_DRMS_M:.2f} m")
    print(f"  Mag nav without TL (est.) : {drms_raw_est:.2f} m")
    print(f"  Mag nav with TL    (est.) : {drms_tl_est:.2f} m")
    improvement = (drms_raw_est - drms_tl_est) / drms_raw_est * 100 \
                  if drms_raw_est > 0 else 0
    print(f"  Estimated improvement     : {improvement:.0f}%")

    if hdg_corr_comp < hdg_corr_raw:
        print("\n  ✓ TL compensation reduced heading correlation — working correctly.")
    else:
        print("\n  ⚠ Heading correlation did not decrease.")
        print("    Possible causes:")
        print("    • Gazebo aircraft model differs significantly from Flt1006")
        print("    • Very short flight (< 2 min) — insufficient heading coverage")
        print("    • Coefficient mismatch — consider re-training on this aircraft")

    # ── Plot ─────────────────────────────────────────────────────────────────
    print(f"\n▶ Generating plot...")
    out_path = make_plot(
        data, Bt_raw, Bt_comp, interference, metrics,
        ulg_name=ulg_path.name, output_dir=out_dir
    )
    print(f"  ✓ Saved: {out_path}")

    print("\n" + "=" * 60)
    print("  To compute exact navigation DRMS:")
    print("    1. Convert this log to HDF5:")
    print(f"       julia gazebo/scripts/ulg_to_xyz.jl {ulg_path}")
    print("    2. Run MagNav.jl benchmark:")
    print("       julia --project=benchmark benchmark/bench_tl_compare.jl")
    print("=" * 60)


if __name__ == "__main__":
    main()
