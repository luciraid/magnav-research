#!/usr/bin/env julia
"""
generate_plots.jl

Generate sample plots demonstrating all visualization functions
in the MagNav TL Benchmark Framework.
"""

using Plots
using Statistics
using Random
using Printf

# Define required structs for demonstration
struct BenchmarkResult
    label::String
    t::Vector{Float64}
    err_north_m::Vector{Float64}
    err_east_m::Vector{Float64}
end

struct NavMetrics
    label       :: String
    rms_north   :: Float64      # m
    rms_east    :: Float64      # m
    drms        :: Float64      # m  √(σN² + σE²)
    drms_2      :: Float64      # m  2×DRMS
    cep         :: Float64      # m  50th percentile of radial error
    final_error :: Float64      # m  at last timestep
    peak_error  :: Float64      # m  maximum radial error
    hdg_corr    :: Float64      # |corr(mag_residual, heading)|
    mag_rms     :: Float64      # nT residual RMS of compensated signal
end

# Include the plotting functions
include("bench_plots.jl")
include("nav_metrics.jl")

# Set up GR backend for better compatibility
gr()

# Create synthetic benchmark results for demonstration
struct BenchmarkResult
    label::String
    t::Vector{Float64}
    err_north_m::Vector{Float64}
    err_east_m::Vector{Float64}
end

# Generate synthetic data
Random.seed!(42)
n_points = 1000
t = range(0, 180, length=n_points)  # 3 minutes at 10Hz

# Create two synthetic results (official vs custom TL)
err_north_official = 0.5 * randn(n_points) .+ 0.1 * sin.(t/10)
err_east_official = 0.3 * randn(n_points) .+ 0.05 * cos.(t/15)

err_north_custom = 0.3 * randn(n_points) .+ 0.05 * sin.(t/10)
err_east_custom = 0.2 * randn(n_points) .+ 0.03 * cos.(t/15)

result_official = BenchmarkResult("Official TL", t, err_north_official, err_east_official)
result_custom = BenchmarkResult("Custom 9-term TL", t, err_north_custom, err_east_custom)

results = [result_official, result_custom]

println("🎨 Generating sample plots for MagNav TL Benchmark Framework...")
println("="^60)

# Create results directory
mkpath("results/plots_demo")

# 1. Error vs Time Plot
println("📊 Generating Error vs Time plot...")
p1 = plot_error_vs_time(results; save_path="results/plots_demo/error_vs_time.png")
println("   ✅ Saved: results/plots_demo/error_vs_time.png")

# 2. Position Error Scatter Plot
println("📊 Generating Position Error Scatter plot...")
p2 = plot_position_errors(results; save_path="results/plots_demo/position_scatter.png")
println("   ✅ Saved: results/plots_demo/position_scatter.png")

# 3. Heading Correlation Plot (simplified version)
println("📊 Generating Heading Correlation plot...")
# Create synthetic heading data
heading = mod.(t * 2, 360)  # Heading changes over time
radial_error_official = sqrt.(err_north_official.^2 .+ err_east_official.^2)
radial_error_custom = sqrt.(err_north_custom.^2 .+ err_east_custom.^2)

# Simple scatter plot instead of full correlation analysis
p3 = plot(heading, radial_error_official;
          title = "Heading vs Radial Error",
          xlabel = "Heading (°)",
          ylabel = "Radial Error (m)",
          label = "Official TL",
          color = COL_OFFICIAL,
          alpha = 0.6,
          PLOT_DEFAULTS...,
          size = (800, 600))

plot!(p3, heading, radial_error_custom;
      label = "Custom 9-term TL",
      color = COL_CUSTOM,
      alpha = 0.6)

savefig(p3, "results/plots_demo/heading_correlation.png")
println("   ✅ Saved: results/plots_demo/heading_correlation.png")

# 4. Magnetic Field Comparison
println("📊 Generating Magnetic Field Comparison plot...")
# Create synthetic magnetic field data
mag_uncomp = 50000 .+ 100 * sin.(t/20) .+ 50 * randn(n_points)
mag_comp_official = mag_uncomp .* 0.7 .+ 30 * randn(n_points)
mag_comp_custom = mag_uncomp .* 0.6 .+ 20 * randn(n_points)

p4 = plot(t, mag_uncomp;
          title = "Magnetic Field Compensation Comparison",
          xlabel = "Time (s)",
          ylabel = "Magnetic Field (nT)",
          label = "Uncompensated",
          color = :gray,
          alpha = 0.7,
          PLOT_DEFAULTS...,
          size = (1000, 600))

plot!(p4, t, mag_comp_official;
      label = "Official TL Compensated",
      color = COL_OFFICIAL,
      lw = 2)

plot!(p4, t, mag_comp_custom;
      label = "Custom 9-term TL Compensated",
      color = COL_CUSTOM,
      lw = 2)

# Add RMS annotations
rms_uncomp = sqrt(mean(mag_uncomp.^2))
rms_official = sqrt(mean(mag_comp_official.^2))
rms_custom = sqrt(mean(mag_comp_custom.^2))

annotate!(p4, 20, maximum(mag_uncomp)*0.95,
         text(@sprintf("Uncomp RMS: %.0f nT", rms_uncomp), FONT_SIZE_ANNOT, :gray, :left))
annotate!(p4, 20, maximum(mag_uncomp)*0.90,
         text(@sprintf("Official RMS: %.0f nT (%.1f%% reduction)", rms_official,
                      100*(1-rms_official/rms_uncomp)), FONT_SIZE_ANNOT, COL_OFFICIAL, :left))
annotate!(p4, 20, maximum(mag_uncomp)*0.85,
         text(@sprintf("Custom RMS: %.0f nT (%.1f%% reduction)", rms_custom,
                      100*(1-rms_custom/rms_uncomp)), FONT_SIZE_ANNOT, COL_CUSTOM, :left))

savefig(p4, "results/plots_demo/magnetic_comparison.png")
println("   ✅ Saved: results/plots_demo/magnetic_comparison.png")

# 5. Residual Analysis Plot (simplified without qqplot)
println("📊 Generating Residual Analysis plot...")
# Calculate statistics for residual analysis
rms_uncomp = sqrt(mean(mag_uncomp.^2))
rms_comp = sqrt(mean(mag_comp_custom.^2))
improvement = (1 - rms_comp/rms_uncomp) * 100
corr_uncomp = cor(heading, mag_uncomp)
corr_comp = cor(heading, mag_comp_custom)

# Create a simplified residual analysis without qqplot
p5_1 = plot(heading, mag_uncomp;
            title = @sprintf("Uncompensated\nRMS: %.1f nT", rms_uncomp),
            xlabel = "Heading (°)",
            ylabel = "Magnetic Residual (nT)",
            label = "Uncompensated",
            color = COL_OFFICIAL,
            PLOT_DEFAULTS...,
            fontsize = FONT_SIZE_LABEL)

plot!(p5_1, heading, mag_comp_custom;
      label = @sprintf("Compensated\nRMS: %.1f nT (%.1f%% improvement)",
                      rms_comp, improvement),
      color = COL_CUSTOM,
      lw = 2)

# Residual histogram
p5_2 = histogram(mag_uncomp;
                 title = "Residual Distribution",
                 xlabel = "Magnetic Residual (nT)",
                 ylabel = "Frequency",
                 label = "Uncompensated",
                 color = COL_OFFICIAL,
                 alpha = 0.7,
                 PLOT_DEFAULTS...,
                 fontsize = FONT_SIZE_LABEL)

histogram!(p5_2, mag_comp_custom;
           label = "Compensated",
           color = COL_CUSTOM,
           alpha = 0.7)

# Simple correlation plot instead of QQ plot
p5_3 = plot(mag_uncomp, mag_comp_custom;
            title = "Residual Correlation",
            xlabel = "Uncompensated (nT)",
            ylabel = "Compensated (nT)",
            label = "Data points",
            color = COL_CUSTOM,
            alpha = 0.6,
            PLOT_DEFAULTS...,
            fontsize = FONT_SIZE_LABEL)

# Add correlation line
corr_val = cor(mag_uncomp, mag_comp_custom)
plot!(p5_3, mag_uncomp, mag_uncomp;
      label = @sprintf("Perfect correlation\n(r = %.3f)", corr_val),
      color = :gray,
      ls = :dash,
      lw = 2)

# Statistics panel
p5_4 = plot(title = "Compensation Statistics",
            grid = false, axis = false, legend = false)

stats_text = """
RMS Improvement: $(round(improvement; digits=1))%
Heading Correlation:
  Uncomp: $(round(corr_uncomp; digits=3))
  Comp: $(round(corr_comp; digits=3))
Variance Reduction: $(round((1 - var(mag_comp_custom)/var(mag_uncomp))*100; digits=1))%
"""

annotate!(p5_4, 0.1, 0.9, text(stats_text, FONT_SIZE_ANNOT, :left))

p5 = plot(p5_1, p5_2, p5_3, p5_4;
          layout = (2,2),
          size = (1200, 900),
          dpi = 300,
          plot_title = "TL Compensation Residual Analysis",
          plot_titlefontsize = FONT_SIZE_TITLE)

savefig(p5, "results/plots_demo/residual_analysis.png")
println("   ✅ Saved: results/plots_demo/residual_analysis.png")

# 6. Benchmark Summary Dashboard
println("📊 Generating Benchmark Summary dashboard...")
p6 = plot(p1, p2, p3, p5;
          layout = (2, 2),
          size = (1400, 1000),
          dpi = 300,
          plot_title = "MagNav TL Benchmark Framework - Results Summary",
          plot_titlefontsize = FONT_SIZE_TITLE + 2)

savefig(p6, "results/plots_demo/benchmark_dashboard.png")
println("   ✅ Saved: results/plots_demo/benchmark_dashboard.png")

println("="^60)
println("🎉 All plots generated successfully!")
println()
println("📁 Plot files saved in: results/plots_demo/")
println()
println("Generated plots:")
println("  • error_vs_time.png        - Position error time series")
println("  • position_scatter.png     - 2D error scatter plot")
println("  • heading_correlation.png  - Heading vs radial error")
println("  • magnetic_comparison.png  - Compensation performance")
println("  • residual_analysis.png    - Statistical analysis")
println("  • benchmark_dashboard.png  - Complete results overview")
println()
println("📊 Key Results Demonstrated:")
println("  • Custom 9-term TL: ~40% RMS reduction")
println("  • Official TL: ~30% RMS reduction")
println("  • Position accuracy: Sub-meter level performance")
println("  • Heading correlation: Improved magnetic stability")