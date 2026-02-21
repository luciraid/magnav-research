#!/usr/bin/env julia
"""
    experiments/bench_multi_segment.jl

Multi-segment cross-validation benchmark.

Runs the Official vs Custom TL comparison across ALL defined SegmentSpecs,
aggregating metrics to produce statistically robust comparison tables and
segment-comparison plots.

This is the primary experiment for publication-quality results.

Usage:
    julia --project=. experiments/bench_multi_segment.jl
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

# -- Config --------------------------------------------------------------------
# Automatically find HDF5 via MagNav's artifact system
const H5_DIR       = sgl_2020_train()  # Downloads artifact on first call
const H5_PATH      = joinpath(H5_DIR, "Flt1006_train.h5")
const MAP_ARTIFACT = "ottawa_area_maps"
const RESULTS_DIR  = get(ENV, "RESULTS_DIR", joinpath(@__DIR__, "..", "results"))

cfg = default_config(; results_dir=RESULTS_DIR, notes="multi_segment_xval")
run_dir = joinpath(RESULTS_DIR, "multiseg_$(cfg.run_id)")
mkpath(joinpath(run_dir, "figures"))
mkpath(joinpath(run_dir, "tables"))
save_config(cfg, joinpath(run_dir, "config.txt"))

println("=" ^ 60)
println("  Multi-Segment Cross-Validation Benchmark")
println("  $(now())")
println("=" ^ 60)

# -- Load data -----------------------------------------------------------------
data   = load_benchmark_data(H5_PATH, MAP_ARTIFACT; verbose=true)
xyz    = data.xyz
mapS   = data.mapS
specs  = sgl2020_flt1006_segments()
segs   = select_segments(xyz, specs; verbose=true)

# -- Accumulate results across segments ----------------------------------------
segment_metric_pairs = NamedTuple[]

for seg in segs
    cal_ind  = seg.cal_ind
    eval_ind = seg.eval_ind
    spec     = seg.spec

    println("\n" * "-" ^ 50)
    @printf "  Segment: %s\n" spec.label

    # Fit custom model on this segment's calibration window
    cal_sig = extract_cal_signals(xyz, cal_ind;
                                   mag_field=cfg.mag_use,
                                   flux_field=cfg.flux_use)
    model = CustomTLModel()
    fit_custom_tl!(model,
                   cal_sig.Bx, cal_sig.By, cal_sig.Bz,
                   cal_sig.mag_scalar, cal_sig.heading;
                   lambda=cfg.tl_lambda, detrend=true, verbose=false)

    @printf "  Custom TL fit: VR=%.1f%%  RMS=%.3f nT  |r_hdg|: %.4f→%.4f\n" \
        100*model.var_reduction model.residual_rms \
        model.heading_corr_pre model.heading_corr_post

    # Run EKF for both methods
    res = run_both_benchmarks(
        data, cal_ind, eval_ind, model, mapS;
        config=cfg, verbose=false
    )

    m_off  = compute_nav_metrics(res.official)
    m_cust = compute_nav_metrics(res.custom)

    @printf "  Official:  DRMS=%.1fm  N=%.1fm  E=%.1fm  |r_hdg|=%.4f\n" \
        m_off.drms m_off.rms_north m_off.rms_east m_off.hdg_corr
    @printf "  Custom:    DRMS=%.1fm  N=%.1fm  E=%.1fm  |r_hdg|=%.4f\n" \
        m_cust.drms m_cust.rms_north m_cust.rms_east m_cust.hdg_corr

    Δ = delta_metrics(m_off, m_cust)
    @printf "  Δ DRMS: %+.1f m (%+.1f%%)\n" Δ.Δdrms Δ.pct_drms

    # Save per-segment figures
    seg_fig_dir = joinpath(run_dir, "figures", spec.label)
    save_all_plots([res.official, res.custom], [m_off, m_cust], seg_fig_dir;
                   tag=spec.label)

    push!(segment_metric_pairs, (
        official = m_off,
        custom   = m_cust,
        label    = spec.label
    ))
end

# -- Aggregate statistics ------------------------------------------------------
println("\n" * "=" ^ 60)
println("  Aggregate Statistics Across $(length(segment_metric_pairs)) Segments")
println("=" ^ 60)

for metric in [:drms, :rms_north, :rms_east, :cep, :final_error, :hdg_corr]
    vals_off  = [getfield(p.official, metric) for p in segment_metric_pairs]
    vals_cust = [getfield(p.custom,   metric) for p in segment_metric_pairs]
    Δvals     = vals_cust .- vals_off

    m_name = rpad(string(metric), 14)
    @printf "  %s  Official: μ=%.2f σ=%.2f  |  Custom: μ=%.2f σ=%.2f  |  Δμ=%+.2f\n" \
        m_name mean(vals_off) std(vals_off) mean(vals_cust) std(vals_cust) mean(Δvals)
end

# -- Aggregate metrics CSV -----------------------------------------------------
rows = []
for p in segment_metric_pairs
    for (method, m) in [("official_tl", p.official), ("custom_tl", p.custom)]
        push!(rows, (
            segment    = p.label,
            method     = method,
            rms_north  = m.rms_north,
            rms_east   = m.rms_east,
            drms       = m.drms,
            cep        = m.cep,
            final_m    = m.final_error,
            peak_m     = m.peak_error,
            hdg_corr   = m.hdg_corr,
            mag_rms_nT = m.mag_rms
        ))
    end
end
agg_df = DataFrame(rows)
agg_path = joinpath(run_dir, "tables", "aggregate_metrics.csv")
CSV.write(agg_path, agg_df)
println("\n  Aggregate CSV: $agg_path")

# -- Multi-segment comparison plots -------------------------------------------
for metric in [:drms, :rms_north, :hdg_corr]
    p = plot_multi_segment(segment_metric_pairs; metric=metric,
            save_path=joinpath(run_dir, "figures", "multiseg_$(metric).png"))
    println("  Figure: multiseg_$metric.png")
end

println("\n" * "=" ^ 60)
@printf "  Multi-segment run complete: %s\n" cfg.run_id
@printf "  Output dir: %s\n" run_dir
println("=" ^ 60)
