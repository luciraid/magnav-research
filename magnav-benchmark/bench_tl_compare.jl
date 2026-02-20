#!/usr/bin/env julia
"""
    experiments/bench_tl_compare.jl

Main benchmark experiment: Official TL vs Custom TL, single-segment.

The HDF5 data file (Flt1006_train.h5) is downloaded automatically via MagNav's
artifact system on first run. No manual data setup required.

Usage (from project root):
    julia --project=. bench_tl_compare.jl

Or with custom output directory:
    RESULTS_DIR=/tmp/results julia --project=. bench_tl_compare.jl

Environment variables
---------------------
RESULTS_DIR : output directory           (default: results/)
RUN_NOTES   : free-text annotation       (optional)

Data Location
-------------
The artifact system automatically downloads to: ~/.julia/artifacts/<hash>/sgl_2020_train/
Do not move or manually manage the HDF5 files - let MagNav handle it.
"""

# -- Environment ---------------------------------------------------------------
# Don't use Pkg.activate with pipe - use global environment where packages installed

push!(LOAD_PATH, @__DIR__)

# MagNav and core libraries
using MagNav
using Printf
using Dates
using LinearAlgebra
using Statistics
using StatsBase
using DataFrames

# Include framework modules directly
include("custom_tl.jl")
include("tl_injection.jl")
include("ekf_wrapper.jl")
include("nav_metrics.jl")
include("bench_plots.jl")
include("segment_utils.jl")
include("data_loader.jl")
include("experiment_config.jl")

# -- Configuration -------------------------------------------------------------
# Automatically find HDF5 via MagNav's artifact system
const H5_DIR      = sgl_2020_train()  # Downloads artifact on first call
const H5_PATH     = joinpath(H5_DIR, "Flt1006_train.h5")
const MAP_ARTIFACT = "ottawa_area_maps"
const RESULTS_DIR  = get(ENV, "RESULTS_DIR", joinpath(@__DIR__, "..", "results"))
const RUN_NOTES    = get(ENV, "RUN_NOTES",   "bench_tl_compare single-segment")

cfg = default_config(; results_dir=RESULTS_DIR, notes=RUN_NOTES)
cfg.detrend_window = :local   # safer for cross-segment evaluation

run_dir = joinpath(RESULTS_DIR, "run_$(cfg.run_id)")
mkpath(joinpath(run_dir, "figures"))
mkpath(joinpath(run_dir, "tables"))
mkpath(joinpath(run_dir, "logs"))

save_config(cfg, joinpath(run_dir, "config.txt"))
describe_config(cfg)

# -- Load data -----------------------------------------------------------------
println("=" ^ 60)
println("  MagNav TL Benchmark: Official vs Custom")
println("  $(now())")
println("=" ^ 60)

data = load_benchmark_data(H5_PATH, MAP_ARTIFACT; verbose=true)
describe_data(data)
xyz  = data.xyz
mapS = data.mapS

# -- Select segments -----------------------------------------------------------
println("\n-- Segment validation -----------------------------------")
specs    = sgl2020_flt1006_segments()
segs     = select_segments(xyz, specs; verbose=true)

# Use first valid segment for the primary run
seg      = segs[1]
cal_ind  = seg.cal_ind
eval_ind = seg.eval_ind
spec     = seg.spec

@printf "\nPrimary segment: %s\n" spec.label
@printf "  Calibration : %d samples\n" sum(cal_ind)
@printf "  Evaluation  : %d samples\n" sum(eval_ind)

# -- Fit Custom TL model -------------------------------------------------------
println("\n-- Fitting Custom TL (9-term) ---------------------------")
cal_signals = extract_cal_signals(xyz, cal_ind;
                                   mag_field  = cfg.mag_use,
                                   flux_field = cfg.flux_use)

custom_model = CustomTLModel()
fit_custom_tl!(custom_model,
               cal_signals.Bx, cal_signals.By, cal_signals.Bz,
               cal_signals.mag_scalar, cal_signals.heading;
               lambda  = cfg.tl_lambda,
               detrend = true,
               verbose = true)

# Save model coefficients
coef_path = joinpath(run_dir, "logs", "custom_tl_coef.txt")
open(coef_path, "w") do io
    println(io, "# Custom 9-term TL coefficients — $(now())")
    println(io, "# Segment: $(spec.label)  |  cal samples: $(custom_model.n_samples)")
    for (i, c) in enumerate(custom_model.coef)
        @printf io "c[%02d] = %.8f\n" i c
    end
    @printf io "var_reduction   = %.4f\n" custom_model.var_reduction
    @printf io "residual_rms    = %.4f nT\n" custom_model.residual_rms
    @printf io "hdg_corr_pre    = %.6f\n" custom_model.heading_corr_pre
    @printf io "hdg_corr_post   = %.6f\n" custom_model.heading_corr_post
end
println("  Coefficients saved: $coef_path")

# -- Run both EKF benchmarks ---------------------------------------------------
println("\n-- Running EKF benchmarks -------------------------------")
results = run_both_benchmarks(
    data, cal_ind, eval_ind, custom_model, mapS;
    config  = cfg,
    verbose = true
)
result_official = results.official
result_custom   = results.custom

all_results = [result_official, result_custom]

# -- Compute metrics -----------------------------------------------------------
println("\n-- Navigation Metrics -----------------------------------")
m_official = compute_nav_metrics(result_official)
m_custom   = compute_nav_metrics(result_custom)
all_metrics = [m_official, m_custom]

print_metrics_table(all_metrics)

# Delta summary
Δ = delta_metrics(m_official, m_custom)
@printf "\n-- Δ Custom − Official ---------------------------------\n"
@printf "  ΔDRMS        : %+.2f m  (%+.1f%%)\n"  Δ.Δdrms     Δ.pct_drms
@printf "  ΔRMS_North   : %+.2f m\n"              Δ.Δrms_north
@printf "  ΔRMS_East    : %+.2f m\n"              Δ.Δrms_east
@printf "  ΔFinal error : %+.2f m\n"              Δ.Δfinal_error
@printf "  Δ|hdg corr|  : %+.5f\n"               Δ.Δhdg_corr
@printf "--------------------------------------------------------\n"

# Save metrics CSV
using CSV, DataFrames
csv_path = joinpath(run_dir, "tables", "metrics_$(spec.label).csv")
CSV.write(csv_path, metrics_table(all_metrics))
println("\n  Metrics saved: $csv_path")

# -- Plots ---------------------------------------------------------------------
println("\n-- Generating figures -----------------------------------")
fig_dir = joinpath(run_dir, "figures")

save_all_plots(all_results, all_metrics, fig_dir; tag=spec.label)

# -- Done ----------------------------------------------------------------------
println("\n" * "=" ^ 60)
@printf "  Run complete: %s\n" cfg.run_id
@printf "  Output dir  : %s\n" run_dir
println("=" ^ 60)
