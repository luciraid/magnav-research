"""
    BenchmarkFramework

Top-level module for the MagNav TL compensation benchmark suite.

Provides:
- Custom 9-term Tolles-Lawson compensation
- Navigation-level side-by-side benchmarking against MagNav official TL
- Metrics: RMS North/East, DRMS, final position error, heading correlation
- Reproducible experiment structure with segment control

Usage:
    julia --project=. src/BenchmarkFramework.jl
    julia --project=. experiments/bench_tl_compare.jl
"""
module BenchmarkFramework

using MagNav
using LinearAlgebra
using Statistics
using StatsBase
using DataFrames
using Dates
using Printf
using Plots

# -- Sub-modules -------------------------------------------------------------
include("custom_tl.jl")
include("tl_injection.jl")
include("ekf_wrapper.jl")
include("nav_metrics.jl")
include("bench_plots.jl")
include("segment_utils.jl")
include("data_loader.jl")
include("experiment_config.jl")

# -- Public API ---------------------------------------------------------------
export
    # TL
    CustomTLModel,
    fit_custom_tl!,
    compensate_custom_tl,
    fit_official_tl,
    compensate_official_tl,
    inject_compensation,

    # EKF wrapper
    run_ekf_benchmark,
    BenchmarkResult,

    # Metrics
    compute_nav_metrics,
    NavMetrics,
    heading_correlation,

    # Plots
    plot_position_errors,
    plot_error_vs_time,
    plot_heading_correlation,
    plot_benchmark_summary,

    # Segments / data
    load_benchmark_data,
    BenchmarkData,
    select_segments,
    SegmentSpec,

    # Config
    ExperimentConfig,
    default_config,
    save_config,
    load_config

end # module
