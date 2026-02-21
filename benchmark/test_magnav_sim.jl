#!/usr/bin/env julia
"""
    benchmark/test_magnav_sim.jl

Run MagNav.jl's built-in synthetic flight simulation and compare results
against the real Flt1006 benchmark results.

Usage:
    cd ~/magnav/magnav-research
    julia benchmark/test_magnav_sim.jl

What this does:
    1. Generates a synthetic 10-minute flight using create_XYZ0() on the
       built-in NAMAD map (no data download required)
    2. Runs MagNav's EKF navigation with their default compensated mag signal
    3. Computes CRLB, INS, and filter DRMS
    4. Prints a side-by-side comparison with the real Flt1006 results

Reference: ~/magnav/MagNav.jl/examples/pluto_sim.jl
"""

push!(LOAD_PATH, @__DIR__)

using MagNav
using Random: seed!
using Statistics: mean, std
using Printf

# ── Real Flt1006 benchmark results (run_20260220_221033) ─────────────────────
# These are the numbers from bench_tl_compare.jl on real SGL 2020 flight data.
const REAL_DRMS_OFFICIAL = 0.59   # m  — MagNav official TL
const REAL_DRMS_CUSTOM   = 0.36   # m  — custom 9-term TL (this project)
const REAL_IMPROVEMENT   = round((1 - REAL_DRMS_CUSTOM / REAL_DRMS_OFFICIAL) * 100; digits=0)

# ── Helper: DRMS formula (matches MagNav.jl convention) ──────────────────────
drms(n_err, e_err) = sqrt(mean(n_err.^2 .+ e_err.^2))

# ── Separator helpers ─────────────────────────────────────────────────────────
banner(s)  = println("\n" * "═" ^ 60 * "\n  " * s * "\n" * "═" ^ 60)
section(s) = println("\n── " * s * " " * "─" ^ max(0, 54 - length(s)))

# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

banner("MagNav Synthetic Flight Simulation")

# ── 1. Load built-in NAMAD map ────────────────────────────────────────────────
section("Loading NAMAD map (built-in, no download)")
mapS_full = get_map(MagNav.namad)
lat_min, lat_max = extrema(rad2deg.(mapS_full.yy))
lon_min, lon_max = extrema(rad2deg.(mapS_full.xx))
@printf "  Map loaded: %.2f°–%.2f° lat, %.2f°–%.2f° lon\n" lat_min lat_max lon_min lon_max

# ── 2. Generate synthetic 10-minute flight ────────────────────────────────────
section("Generating synthetic 10-minute flight")
seed!(33)                                        # matches pluto_sim.jl for reproducibility
t   = 600                                        # seconds
xyz = create_XYZ0(mapS_full; alt=mapS_full.alt, t=t)
traj = xyz.traj
ins  = xyz.ins

n_samples = length(traj.tt)
@printf "  Samples : %d  (%.1f Hz)\n"            n_samples (1/traj.dt)
@printf "  Duration: %.0f s  (%.1f min)\n"       t (t/60)
@printf "  Altitude  : %.0f m\n"                   traj.alt[1]
@printf "  Lat range : %.4f°–%.4f°\n"             rad2deg(minimum(traj.lat)) rad2deg(maximum(traj.lat))
@printf "  Lon range : %.4f°–%.4f°\n"             rad2deg(minimum(traj.lon)) rad2deg(maximum(traj.lon))

# ── 3. Trim map and build interpolant ─────────────────────────────────────────
section("Preparing map interpolant")
mapS     = map_trim(mapS_full, traj; pad=10)
itp_mapS = map_interpolate(mapS)
@printf "  Trimmed map: %d × %d grid points\n"  size(mapS.map, 1) size(mapS.map, 2)

# ── 4. Build EKF model parameters ────────────────────────────────────────────
section("Creating EKF model")
(P0, Qd, R) = create_model(traj.dt, traj.lat[1])
@printf "  Measurement noise R = %.1f nT²\n"    R

# ── 5. Run EKF navigation ─────────────────────────────────────────────────────
section("Running EKF navigation (MagNav default TL, mag_1_c)")
println("  This may take 30–60 s on first run (NAMAD map interpolation)...")
flush(stdout)

mag_use = xyz.mag_1_c   # MagNav's built-in compensated scalar magnetometer

t_start = time()
(crlb_out, ins_out, filt_out) = run_filt(traj, ins, mag_use, itp_mapS, :ekf; P0, Qd, R)
elapsed = time() - t_start
@printf "  EKF completed in %.1f s\n"           elapsed

# ── 6. Compute DRMS metrics ───────────────────────────────────────────────────
section("Computing DRMS metrics")

crlb_drms = drms(crlb_out.n_std, crlb_out.e_std)
ins_drms  = drms(ins_out.n_err,  ins_out.e_err)
filt_drms = drms(filt_out.n_err, filt_out.e_err)

@printf "  CRLB (theoretical min) : %7.3f m\n"  crlb_drms
@printf "  INS  (dead reckoning)  : %7.3f m\n"  ins_drms
@printf "  EKF filter             : %7.3f m\n"  filt_drms

# ── 7. Comparison table ───────────────────────────────────────────────────────
banner("Results Comparison")

println("""
┌──────────────────────────────────────────────┬──────────────┐
│ Method                                       │  DRMS        │
├──────────────────────────────────────────────┼──────────────┤""")

# Synthetic flight rows
@printf "│  Synthetic (NAMAD) — CRLB lower bound        │  %8.3f m  │\n"  crlb_drms
@printf "│  Synthetic (NAMAD) — INS dead reckoning      │  %8.1f m  │\n"   ins_drms
@printf "│  Synthetic (NAMAD) — EKF (MagNav default TL) │  %8.3f m  │\n"  filt_drms

println("├──────────────────────────────────────────────┼──────────────┤")

# Real Flt1006 rows
@printf "│  Real Flt1006      — MagNav official TL      │  %8.2f m  │\n"  REAL_DRMS_OFFICIAL
@printf "│  Real Flt1006      — Custom 9-term TL ★      │  %8.2f m  │\n"  REAL_DRMS_CUSTOM

println("└──────────────────────────────────────────────┴──────────────┘")
println("  ★  This project (run_20260220_221033)")

# ── 8. Interpretation ─────────────────────────────────────────────────────────
section("Interpretation")

sim_vs_real = filt_drms / REAL_DRMS_OFFICIAL
custom_vs_official = REAL_DRMS_CUSTOM / REAL_DRMS_OFFICIAL

println()
@printf "  Synthetic NAMAD baseline  :  %.3f m DRMS\n"         filt_drms
@printf "  Real Flt1006 official TL  :  %.2f m DRMS\n"          REAL_DRMS_OFFICIAL
@printf "  Real Flt1006 custom TL ★  :  %.2f m DRMS  (%.0f%% better than official)\n" REAL_DRMS_CUSTOM REAL_IMPROVEMENT

println()

if filt_drms < REAL_DRMS_CUSTOM
    @printf "  ✓ Synthetic NAMAD result (%.3f m) is better than our real Flt1006\n"  filt_drms
    println("    custom TL result — expected: NAMAD is a smooth, high-quality map")
    println("    purpose-built for simulation. Real SGL data has map mismatch,")
    println("    sensor noise, and genuine aircraft interference.")
elseif filt_drms < REAL_DRMS_OFFICIAL
    @printf "  ✓ Synthetic NAMAD result (%.3f m) falls between our custom TL\n"  filt_drms
    println("    and official TL on real flight data — consistent with expectations.")
else
    @printf "  ℹ Synthetic NAMAD result (%.3f m) is higher than real Flt1006\n"  filt_drms
    println("    results. This can occur if the NAMAD map has low gradient in the")
    println("    generated flight region (less information for the EKF).")
end

println()
println("  Key takeaway:")
println("  The NAMAD CRLB ($(round(Int,crlb_drms)) m) and Flt1006 results (0.36 m) are not directly")
println("  comparable — NAMAD is a coarse continental map; the Ottawa SGL maps used")
println("  in Flt1006 have much higher spatial resolution and gradient, enabling")
println("  sub-metre navigation. The CRLB for the Ottawa map would be far lower.")
@printf "  EKF on NAMAD:            %6.2f m DRMS  (constrained by map quality)\n" filt_drms
@printf "  Our custom TL, Flt1006:  %6.2f m DRMS  (high-res Ottawa SGL map)\n"   REAL_DRMS_CUSTOM

println()
println("  To run the full real-data benchmark:")
println("    julia --project=benchmark benchmark/bench_tl_compare.jl")
println()
