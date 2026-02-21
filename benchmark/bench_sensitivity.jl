#!/usr/bin/env julia
"""
    experiments/bench_sensitivity.jl

Sensitivity analysis: how does custom TL performance vary with calibration
window length and heading diversity?

This experiment answers: "How much calibration data is enough?"

Sweeps calibration window duration from 60s to 900s using a fixed evaluation
segment, and reports how DRMS and heading correlation change.

Usage:
    julia --project=. experiments/bench_sensitivity.jl
"""

# Don't use Pkg.activate with pipe - use global environment where packages installed

push!(LOAD_PATH, @__DIR__)

using MagNav
using Printf, Dates
using CSV, DataFrames
using LinearAlgebra
using Statistics
using StatsBase

# Include framework modules directly
include("custom_tl.jl")
include("tl_injection.jl")
include("ekf_wrapper.jl")
include("nav_metrics.jl")
include("bench_plots.jl")
include("segment_utils.jl")
include("data_loader.jl")
include("experiment_config.jl")

# Automatically find HDF5 via MagNav's artifact system
const H5_DIR       = sgl_2020_train()  # Downloads artifact on first call
const H5_PATH      = joinpath(H5_DIR, "Flt1006_train.h5")
const MAP_ARTIFACT = "ottawa_area_maps"
const RESULTS_DIR  = get(ENV, "RESULTS_DIR", joinpath(@__DIR__, "..", "results"))

cfg     = default_config(; results_dir=RESULTS_DIR, notes="sensitivity_cal_length")
run_dir = joinpath(RESULTS_DIR, "sensitivity_$(cfg.run_id)")
mkpath(joinpath(run_dir, "figures"))
mkpath(joinpath(run_dir, "tables"))
save_config(cfg, joinpath(run_dir, "config.txt"))

println("=" ^ 60)
println("  Sensitivity Analysis: Calibration Window Length")
println("  $(now())")
println("=" ^ 60)

data = load_benchmark_data(H5_PATH, MAP_ARTIFACT; verbose=true)
xyz  = data.xyz
mapS = data.mapS

# Fixed evaluation segment: [800s, 1400s]
EVAL_START = 800.0
EVAL_END   = 1400.0

# Calibration window starts at 0, sweeps end time
CAL_START     = 0.0
CAL_END_SWEEP = [60.0, 120.0, 180.0, 300.0, 450.0, 600.0, 750.0, 900.0]
GAP           = 60.0   # gap before eval always maintained

t_vec = Float64.(xyz.traj.tt) .- Float64(xyz.traj.tt[1])

rows = []
println("\nSweep:")
for cal_end in CAL_END_SWEEP
    if cal_end + GAP > EVAL_START
        @warn "cal_end=$cal_end + gap=$GAP would overlap eval start=$EVAL_START, skipping"
        continue
    end

    cal_ind  = (t_vec .>= CAL_START)  .& (t_vec .<= cal_end)
    eval_ind = (t_vec .>= EVAL_START) .& (t_vec .<= EVAL_END)

    @printf "  cal=[%.0f–%.0fs] (%d pts)  eval=[%.0f–%.0fs] (%d pts)  " \
        CAL_START cal_end sum(cal_ind) EVAL_START EVAL_END sum(eval_ind)

    # Heading diversity in calibration window
    hdg_range = range_heading(xyz.ins_yaw[cal_ind])

    # Fit custom TL
    cal_sig = extract_cal_signals(xyz, cal_ind;
                                   mag_field=cfg.mag_use, flux_field=cfg.flux_use)
    model = CustomTLModel()
    fit_custom_tl!(model, cal_sig.Bx, cal_sig.By, cal_sig.Bz,
                   cal_sig.mag_scalar, cal_sig.heading;
                   lambda=cfg.tl_lambda, detrend=true, verbose=false)

    # Run EKF (custom only for speed; include official for reference)
    res = run_both_benchmarks(data, cal_ind, eval_ind, model, mapS;
                              config=cfg, verbose=false)
    m_off  = compute_nav_metrics(res.official)
    m_cust = compute_nav_metrics(res.custom)

    @printf "DRMS: off=%.1fm  cust=%.1fm  VR=%.1f%%  hdg_range=%.0f°\n" \
        m_off.drms m_cust.drms 100*model.var_reduction rad2deg(hdg_range)

    push!(rows, (
        cal_end_s    = cal_end,
        cal_n_pts    = sum(cal_ind),
        hdg_range_deg= rad2deg(hdg_range),
        var_reduction= model.var_reduction,
        drms_off     = m_off.drms,
        drms_cust    = m_cust.drms,
        hdg_corr_off = m_off.hdg_corr,
        hdg_corr_cust= m_cust.hdg_corr,
    ))
end

df = DataFrame(rows)
csv_path = joinpath(run_dir, "tables", "sensitivity_cal_length.csv")
CSV.write(csv_path, df)
println("\n  CSV saved: $csv_path")

# -- Plots ---------------------------------------------------------------------
using Plots

p1 = plot(df.cal_end_s, [df.drms_off df.drms_cust];
          label=["Official TL" "Custom TL"],
          xlabel="Calibration window (s)",
          ylabel="DRMS (m)",
          title="DRMS vs Calibration Length",
          marker=:circle, lw=2,
          color=[COL_OFFICIAL COL_CUSTOM],
          grid=true, size=(700,400), dpi=150)
savefig(p1, joinpath(run_dir, "figures", "sensitivity_drms.png"))

p2 = plot(df.cal_end_s, df.hdg_range_deg;
          xlabel="Calibration window (s)",
          ylabel="Heading range (°)",
          title="Heading Diversity vs Calibration Length",
          color=:teal, lw=2, marker=:square,
          grid=true, size=(700,350), dpi=150, legend=false)
savefig(p2, joinpath(run_dir, "figures", "sensitivity_hdg_range.png"))

println("  Figures saved.")

println("\n" * "=" ^ 60)
@printf "  Sensitivity run complete: %s\n" cfg.run_id
println("=" ^ 60)

# -- Helper --------------------------------------------------------------------
function range_heading(yaw_rad)
    hdg = rad2deg.(Float64.(yaw_rad))
    return deg2rad(maximum(hdg) - minimum(hdg))
end
