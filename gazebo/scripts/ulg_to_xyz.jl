#!/usr/bin/env julia
"""
    ulg_to_xyz.jl  —  PX4 CSV → MagNav HDF5 converter

Converts a CSV produced by extract_px4_log.py into an HDF5 file that
MagNav's `h5_to_df` / `get_XYZ` pipeline can read directly.

Pipeline position:
    PX4 Gazebo → extract_px4_log.py → <flight>.csv → ulg_to_xyz.jl → <flight>.h5
                                                                            ↓
                                                            benchmark/bench_tl_compare.jl

Usage:
    cd ~/magnav/magnav-research
    julia gazebo/scripts/ulg_to_xyz.jl gazebo/analysis/flight01.csv

Output:
    gazebo/synthetic_h5/flight01.h5

CSV columns expected (from extract_px4_log.py):
    timestamp   — microseconds since boot
    lat, lon    — degrees
    alt         — meters (above sea level)
    accel_x/y/z — m/s² (body frame)
    gyro_x/y/z  — rad/s (body frame)
    mag_x/y/z   — Gauss (body frame, from vehicle_magnetometer)
    roll, pitch, yaw — radians

MagNav HDF5 fields written:
    tt          — time vector (seconds, 0-based)
    lat, lon    — radians (GPS truth)
    alt         — meters
    vn, ve, vd  — velocity NED (m/s, from finite-diff GPS)
    ins_lat, ins_lon, ins_alt  — same as GPS (Gazebo has no real INS drift)
    ins_vn, ins_ve, ins_vd     — same as GPS velocity
    ins_roll, ins_pitch, ins_yaw — from PX4 attitude estimator (radians)
    flux_a_x/y/z — fluxgate vector (nT), converted from Gauss
    mag_1_uc    — uncompensated scalar (nT) = ‖flux_a‖
    mag_1_c     — placeholder (= mag_1_uc, apply TL in benchmark)

Notes:
  • All data is resampled to TARGET_HZ (default 10 Hz) via linear interpolation
    to match the SGL 2020 dataset cadence.
  • No magnetic anomaly is added here. The benchmark framework will overlay
    the Ottawa map anomaly via MagNav's get_map_val.
  • INS fields are copies of GPS truth because Gazebo SITL does not export
    INS integration residuals. This means navigation error is purely from
    magnetic compensation quality, not INS drift — appropriate for TL testing.
"""

using CSV
using DataFrames
using HDF5
using LinearAlgebra
using Statistics

# ── Configuration ─────────────────────────────────────────────────────────────

const TARGET_HZ  = 10          # Hz — must match SGL dataset
const GAUSS_TO_NT = 1e5        # 1 Gauss = 100,000 nT
const DEG_TO_RAD  = π / 180.0

# ── Helpers ───────────────────────────────────────────────────────────────────

"""
    resample_to_hz(t_src, y_src, hz) -> (t_new, y_new)

Linearly interpolate a signal onto a uniform grid at `hz` Hz.
"""
function resample_to_hz(t_src::Vector{Float64}, y_src::Vector{Float64}, hz::Int)
    t0   = t_src[1]
    t1   = t_src[end]
    dt   = 1.0 / hz
    t_new = collect(t0:dt:t1)
    y_new = Vector{Float64}(undef, length(t_new))
    for (i, t) in enumerate(t_new)
        j = searchsortedlast(t_src, t)
        j = clamp(j, 1, length(t_src) - 1)
        frac = (t - t_src[j]) / (t_src[j+1] - t_src[j])
        y_new[i] = y_src[j] * (1 - frac) + y_src[j+1] * frac
    end
    return t_new, y_new
end

"""
    finite_diff_velocity(lat_rad, lon_rad, alt_m, dt) -> (vn, ve, vd)

Estimate NED velocity from GPS position using central differences.
Units: m/s.
"""
function finite_diff_velocity(lat::Vector{Float64},
                               lon::Vector{Float64},
                               alt::Vector{Float64},
                               dt::Float64)
    n   = length(lat)
    R_E = 6_378_137.0  # Earth equatorial radius (m)

    vn  = zeros(n)
    ve  = zeros(n)
    vd  = zeros(n)

    for i in 2:n-1
        dlat = lat[i+1] - lat[i-1]
        dlon = lon[i+1] - lon[i-1]
        dalt = alt[i+1] - alt[i-1]
        vn[i] =  dlat * R_E                    / (2dt)
        ve[i] =  dlon * R_E * cos(lat[i])      / (2dt)
        vd[i] = -dalt                           / (2dt)
    end
    # Forward/backward diff at endpoints
    vn[1]   = vn[2];   ve[1]   = ve[2];   vd[1]   = vd[2]
    vn[end] = vn[end-1]; ve[end] = ve[end-1]; vd[end] = vd[end-1]

    return vn, ve, vd
end

# ── Main ──────────────────────────────────────────────────────────────────────

function main(args)
    if length(args) < 1
        println("""
        Usage: julia ulg_to_xyz.jl <input.csv> [output.h5]

        Defaults:
          output = gazebo/synthetic_h5/<stem>.h5
        """)
        return
    end

    csv_path = args[1]
    stem     = splitext(basename(csv_path))[1]
    out_dir  = joinpath(@__DIR__, "..", "synthetic_h5")
    mkpath(out_dir)
    h5_path  = length(args) >= 2 ? args[2] : joinpath(out_dir, stem * ".h5")

    println("=" ^ 60)
    println("  ulg_to_xyz — PX4 CSV → MagNav HDF5")
    println("=" ^ 60)
    println("\n▶ Reading CSV: $csv_path")

    df = CSV.read(csv_path, DataFrame)
    println("  Raw rows: $(nrow(df))  columns: $(ncol(df))")

    # ── Convert time to seconds ────────────────────────────────────────────────
    t_us  = Float64.(df.timestamp)
    t_s   = (t_us .- t_us[1]) ./ 1e6   # μs → s, zero-based

    # ── Resample everything to TARGET_HZ ──────────────────────────────────────
    println("\n▶ Resampling to $(TARGET_HZ) Hz...")

    function rsmp(col::Symbol)
        _, y = resample_to_hz(t_s, Float64.(df[!, col]), TARGET_HZ)
        return y
    end

    t_ref, _ = resample_to_hz(t_s, t_s, TARGET_HZ)   # uniform time grid (s)
    n        = length(t_ref)
    dt       = 1.0 / TARGET_HZ

    lat_deg = rsmp(:lat);   lon_deg = rsmp(:lon);   alt_m = rsmp(:alt)
    lat_rad = lat_deg .* DEG_TO_RAD
    lon_rad = lon_deg .* DEG_TO_RAD

    roll  = rsmp(:roll)
    pitch = rsmp(:pitch)
    yaw   = rsmp(:yaw)

    # Magnetometer: Gauss → nT
    Bx_nT = rsmp(:mag_x) .* GAUSS_TO_NT
    By_nT = rsmp(:mag_y) .* GAUSS_TO_NT
    Bz_nT = rsmp(:mag_z) .* GAUSS_TO_NT
    Bt_nT = sqrt.(Bx_nT.^2 .+ By_nT.^2 .+ Bz_nT.^2)

    vn, ve, vd = finite_diff_velocity(lat_rad, lon_rad, alt_m, dt)

    println("  Output samples: $n  (~$(round(n/TARGET_HZ/60; digits=1)) min)")
    println("  Lat range : $(round(minimum(lat_deg); digits=4))° – $(round(maximum(lat_deg); digits=4))°")
    println("  Alt range : $(round(minimum(alt_m); digits=0)) – $(round(maximum(alt_m); digits=0)) m")
    println("  |B| range : $(round(minimum(Bt_nT); digits=0)) – $(round(maximum(Bt_nT); digits=0)) nT")

    # ── Write HDF5 ────────────────────────────────────────────────────────────
    println("\n▶ Writing HDF5: $h5_path")

    h5open(h5_path, "w") do h5
        # Truth / GPS trajectory
        h5["tt"]    = t_ref
        h5["lat"]   = lat_rad
        h5["lon"]   = lon_rad
        h5["alt"]   = alt_m
        h5["vn"]    = vn
        h5["ve"]    = ve
        h5["vd"]    = vd

        # INS (copy of GPS — no real INS in SITL)
        h5["ins_lat"]   = lat_rad
        h5["ins_lon"]   = lon_rad
        h5["ins_alt"]   = alt_m
        h5["ins_vn"]    = vn
        h5["ins_ve"]    = ve
        h5["ins_vd"]    = vd
        h5["ins_roll"]  = roll
        h5["ins_pitch"] = pitch
        h5["ins_yaw"]   = yaw

        # Fluxgate vector sensor (body frame, nT)
        h5["flux_a_x"] = Bx_nT
        h5["flux_a_y"] = By_nT
        h5["flux_a_z"] = Bz_nT

        # Scalar magnetometer (uncompensated = raw magnitude)
        h5["mag_1_uc"] = Bt_nT
        # Compensated placeholder — benchmark will replace this via TL injection
        h5["mag_1_c"]  = Bt_nT

        # Metadata
        attrs(h5)["source"]      = csv_path
        attrs(h5)["target_hz"]   = TARGET_HZ
        attrs(h5)["n_samples"]   = n
        attrs(h5)["created_by"]  = "ulg_to_xyz.jl"
    end

    println("  ✓ Wrote $n samples")
    println("\n  Next steps:")
    println("    1. Run the benchmark on this file:")
    println("       FLIGHT_H5=$h5_path julia benchmark/bench_tl_compare.jl")
    println("    2. Or use validate.jl with h5_file pointing to $h5_path")
    println("\n" * "=" ^ 60)
end

main(ARGS)
