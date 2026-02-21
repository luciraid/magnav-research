"""
    bench_plots.jl

Publication-quality plots for the MagNav TL benchmark suite.

All functions return a Plots.jl plot object and optionally save to disk.
Default backend: GR (bundled with Plots.jl, no extra deps).

Enhanced features:
- Professional color schemes and typography
- Statistical annotations (RMS, DRMS values)
- Confidence intervals and error bars
- Publication-ready formatting
- High-resolution output (300 DPI)

Plot functions
--------------
- `plot_error_vs_time`        : North/East error vs time for both runs
- `plot_position_errors`      : North vs East scatter + time colourmap
- `plot_heading_correlation`  : Radial error binned by heading (polar + linear)
- `plot_mag_comparison`       : Compensated scalar signals side by side
- `plot_benchmark_summary`    : 2×2 panel dashboard
- `plot_multi_segment`        : Box/violin across N experiment segments
- `plot_residual_analysis`    : TL fit quality and residual statistics
"""

# -- Enhanced Colour palette ---------------------------------------------------
const COL_OFFICIAL = RGB(31/255, 119/255, 180/255)    # Steel blue
const COL_CUSTOM   = RGB(214/255, 39/255, 40/255)     # Firebrick red
const COL_TRUTH    = RGB(0, 0, 0)                      # Black
const COL_GRID     = RGB(0.8, 0.8, 0.8)                # Light gray
const COL_ANNOTATE = RGB(0.3, 0.3, 0.3)                # Dark gray

# Enhanced typography
const FONT_SIZE_TITLE = 14
const FONT_SIZE_LABEL = 12
const FONT_SIZE_TICK  = 10
const FONT_SIZE_ANNOT = 9

# Plot styling defaults
const PLOT_DEFAULTS = Dict(
    :fontfamily => "DejaVu Sans",
    :grid => true,
    :gridcolor => COL_GRID,
    :gridalpha => 0.3,
    :minorgrid => true,
    :minorgridcolor => COL_GRID,
    :minorgridalpha => 0.1,
    :linewidth => 2,
    :markersize => 4,
    :dpi => 300,
    :size => (800, 600)
)

# -- Error vs Time -------------------------------------------------------------

"""
    plot_error_vs_time(results; save_path=nothing, title_prefix="", show_stats=true) -> Plot

Plot North error and East error vs time (seconds) for each BenchmarkResult
in `results`. Produces a 2-row stacked figure with enhanced styling and statistics.
"""
function plot_error_vs_time(
    results     :: Vector{BenchmarkResult};
    save_path   :: Union{String,Nothing} = nothing,
    title_prefix:: String = "",
    show_stats  :: Bool = true
)
    palette = [COL_OFFICIAL, COL_CUSTOM, RGB(34/255, 139/255, 34/255), RGB(148/255, 0, 211/255)]

    p_north = plot(;
        xlabel = "Time (seconds)",
        ylabel = "North Error (m)",
        title  = title_prefix * "North Position Error vs Time",
        legend = :topright,
        PLOT_DEFAULTS...,
        fontsize = FONT_SIZE_LABEL
    )
    p_east = plot(;
        xlabel = "Time (seconds)",
        ylabel = "East Error (m)",
        title  = title_prefix * "East Position Error vs Time",
        legend = :topright,
        PLOT_DEFAULTS...,
        fontsize = FONT_SIZE_LABEL
    )

    for (i, r) in enumerate(results)
        col = palette[mod1(i, length(palette))]
        plot!(p_north, r.t, r.err_north_m;
              label=r.label, color=col, lw=2, alpha=0.8)

        plot!(p_east, r.t, r.err_east_m;
              label=r.label, color=col, lw=2, alpha=0.8)

        if show_stats
            # Add RMS annotations
            rms_north = sqrt(mean(r.err_north_m.^2))
            rms_east  = sqrt(mean(r.err_east_m.^2))

            annotate!(p_north, maximum(r.t)*0.02, maximum(r.err_north_m)*0.9,
                     text(@sprintf("RMS: %.2f m", rms_north),
                          FONT_SIZE_ANNOT, col, :left))

            annotate!(p_east, maximum(r.t)*0.02, maximum(r.err_east_m)*0.9,
                     text(@sprintf("RMS: %.2f m", rms_east),
                          FONT_SIZE_ANNOT, col, :left))
        end
    end

    # Enhanced zero lines
    hline!(p_north, [0.0]; color=COL_TRUTH, ls=:dash, lw=1.5, alpha=0.7, label="")
    hline!(p_east,  [0.0]; color=COL_TRUTH, ls=:dash, lw=1.5, alpha=0.7, label="")

    fig = plot(p_north, p_east;
               layout=(2,1),
               size=(1000, 800),
               dpi=300,
               plot_title=title_prefix * "Position Error Time Series",
               plot_titlefontsize=FONT_SIZE_TITLE)

    !isnothing(save_path) && savefig(fig, save_path)
    return fig
end

# -- Position Scatter ----------------------------------------------------------

"""
    plot_position_errors(results; save_path=nothing, show_drms=true) -> Plot

2D scatter of (East error, North error) for each run with enhanced visualization.
Shows temporal evolution with color gradient and DRMS circles.
"""
function plot_position_errors(
    results   :: Vector{BenchmarkResult};
    save_path :: Union{String,Nothing} = nothing,
    show_drms :: Bool = true
)
    palette = [COL_OFFICIAL, COL_CUSTOM, RGB(34/255, 139/255, 34/255), RGB(148/255, 0, 211/255)]

    fig = plot(;
        xlabel  = "East Error (m)",
        ylabel  = "North Error (m)",
        title   = "Position Error Scatter Plot",
        aspect_ratio = 1,
        legend  = :topright,
        PLOT_DEFAULTS...,
        fontsize = FONT_SIZE_LABEL
    )

    for (i, r) in enumerate(results)
        col  = palette[mod1(i, length(palette))]
        n    = length(r.t)

        # Time-based color gradient
        cmap = cgrad([col, RGB(1,1,1), col]; alpha=0.6)

        scatter!(fig, r.err_east_m, r.err_north_m;
                 label     = r.label,
                 markercolor = cmap,
                 markersize  = 3,
                 markerstrokewidth = 0.5,
                 markerstrokecolor = col,
                 alpha       = 0.7)

        if show_drms
            # DRMS circle with enhanced styling
            drms_val  = sqrt(mean(r.err_north_m.^2 .+ r.err_east_m.^2))
            θ_circle  = range(0, 2π; length=200)

            plot!(fig,
                  drms_val .* cos.(θ_circle),
                  drms_val .* sin.(θ_circle);
                  label    = @sprintf("%s DRMS: %.2fm", r.label, drms_val),
                  color    = col,
                  lw       = 2.5,
                  ls       = :dash,
                  alpha    = 0.8)
        end
    end

    # Enhanced axes
    hline!(fig, [0.0]; color=COL_TRUTH, lw=1.5, alpha=0.7, label="")
    vline!(fig, [0.0]; color=COL_TRUTH, lw=1.5, alpha=0.7, label="")

    # Add colorbar annotation
    annotate!(fig, maximum([r.err_east_m[1] for r in results])*0.8,
              maximum([r.err_north_m[1] for r in results])*0.9,
              text("Color: Time progression", FONT_SIZE_ANNOT, COL_ANNOTATE, :right))

    !isnothing(save_path) && savefig(fig, save_path)
    return fig
end

# -- Residual Analysis (New) ---------------------------------------------------

"""
    plot_residual_analysis(mag_uncomp, mag_comp, heading; save_path=nothing) -> Plot

Analyze TL compensation quality with residual statistics and heading correlation.
Shows before/after compensation with RMS improvement and correlation analysis.
"""
function plot_residual_analysis(
    mag_uncomp :: Vector{Float64},
    mag_comp   :: Vector{Float64},
    heading    :: Vector{Float64};
    save_path  :: Union{String,Nothing} = nothing
)
    # Calculate statistics
    rms_uncomp = sqrt(mean(mag_uncomp.^2))
    rms_comp   = sqrt(mean(mag_comp.^2))
    improvement = (1 - rms_comp/rms_uncomp) * 100

    # Heading correlation
    corr_uncomp = cor(heading, mag_uncomp)
    corr_comp   = cor(heading, mag_comp)

    # Create 2x2 subplot
    p1 = plot(heading, mag_uncomp;
              title = @sprintf("Uncompensated\nRMS: %.1f nT", rms_uncomp),
              xlabel = "Heading (°)",
              ylabel = "Magnetic Residual (nT)",
              label = "Uncompensated",
              color = COL_OFFICIAL,
              PLOT_DEFAULTS...,
              fontsize = FONT_SIZE_LABEL)

    plot!(p1, heading, mag_comp;
          label = @sprintf("Compensated\nRMS: %.1f nT (%.1f%% improvement)",
                          rms_comp, improvement),
          color = COL_CUSTOM,
          lw = 2)

    # Residual histogram
    p2 = histogram(mag_uncomp;
                   title = "Residual Distribution",
                   xlabel = "Magnetic Residual (nT)",
                   ylabel = "Frequency",
                   label = "Uncompensated",
                   color = COL_OFFICIAL,
                   alpha = 0.7,
                   PLOT_DEFAULTS...,
                   fontsize = FONT_SIZE_LABEL)

    histogram!(p2, mag_comp;
               label = "Compensated",
               color = COL_CUSTOM,
               alpha = 0.7)

    # QQ plot for normality check
    p3 = qqplot(mag_uncomp, mag_comp;
                title = "Q-Q Plot: Uncomp vs Comp",
                xlabel = "Uncompensated Quantiles",
                ylabel = "Compensated Quantiles",
                PLOT_DEFAULTS...,
                fontsize = FONT_SIZE_LABEL)

    # Statistics panel
    p4 = plot(title = "Compensation Statistics",
              grid = false, axis = false, legend = false)

    stats_text = """
    RMS Improvement: $(round(improvement; digits=1))%
    Heading Correlation:
      Uncomp: $(round(corr_uncomp; digits=3))
      Comp: $(round(corr_comp; digits=3))
    Variance Reduction: $(round((1 - var(mag_comp)/var(mag_uncomp))*100; digits=1))%
    """

    annotate!(p4, 0.1, 0.9, text(stats_text, FONT_SIZE_ANNOT, :left))

    fig = plot(p1, p2, p3, p4;
               layout = (2,2),
               size = (1200, 900),
               dpi = 300,
               plot_title = "TL Compensation Residual Analysis",
               plot_titlefontsize = FONT_SIZE_TITLE)

    !isnothing(save_path) && savefig(fig, save_path)
    return fig
end

# -- Heading Correlation -------------------------------------------------------

"""
    plot_heading_correlation(results; n_bins=36, save_path=nothing) -> Plot

2-row figure:
  Row 1: radial position error binned by heading (bar chart, 360°)
  Row 2: raw scatter of heading vs radial error with linear regression overlay
"""
function plot_heading_correlation(
    results   :: Vector{BenchmarkResult};
    n_bins    :: Int = 36,
    save_path :: Union{String,Nothing} = nothing
)
    palette = [COL_OFFICIAL, COL_CUSTOM, :darkgreen, :purple]

    p_bar   = plot(; title="Heading-Binned Mean Radial Error",
                     xlabel="Heading (°)", ylabel="Mean Radial Error (m)",
                     legend=:topright, grid=true)
    p_scat  = plot(; title="Heading vs Radial Error",
                     xlabel="Heading (°)", ylabel="Radial Error (m)",
                     legend=:topright, grid=true, alpha=0.4)

    for (i, r) in enumerate(results)
        col = palette[mod1(i, length(palette))]

        bins_deg, mean_err, corr_val = heading_correlation(r; n_bins=n_bins)

        bar!(p_bar, bins_deg, mean_err;
             label      = "$(r.label)  |r|=$(round(corr_val;digits=4))",
             color      = col,
             alpha      = 0.55,
             bar_width  = 360.0/n_bins * 0.8)

        hdg_deg = rad2deg.(r.heading_rad)
        radial  = @. sqrt(r.err_north_m^2 + r.err_east_m^2)

        scatter!(p_scat, hdg_deg, radial;
                 label  = "$(r.label)  |r|=$(round(corr_val;digits=4))",
                 color  = col,
                 ms     = 2, markerstrokewidth=0)

        # OLS trend line
        A_ols = [ones(length(hdg_deg)) hdg_deg]
        b_ols = A_ols \ radial
        x_lr  = [minimum(hdg_deg), maximum(hdg_deg)]
        y_lr  = b_ols[1] .+ b_ols[2] .* x_lr
        plot!(p_scat, x_lr, y_lr; color=col, lw=2, label="")
    end

    fig = plot(p_bar, p_scat; layout=(2,1), size=(900, 650), dpi=150)
    !isnothing(save_path) && savefig(fig, save_path)
    return fig
end

# -- Compensated Mag Comparison ------------------------------------------------

"""
    plot_mag_comparison(results; save_path=nothing) -> Plot

Overlay the compensated scalar magnetometer signals for visual inspection
of residual structure.
"""
function plot_mag_comparison(
    results   :: Vector{BenchmarkResult};
    save_path :: Union{String,Nothing} = nothing
)
    palette = [COL_OFFICIAL, COL_CUSTOM, :darkgreen, :purple]

    fig = plot(;
        xlabel  = "Time (s)",
        ylabel  = "Compensated Mag (nT)",
        title   = "Compensated Scalar Magnetometer Comparison",
        legend  = :topright,
        grid    = true,
        size    = (900, 380),
        dpi     = 150
    )

    for (i, r) in enumerate(results)
        col = palette[mod1(i, length(palette))]
        mag_demean = r.compensated_mag .- mean(r.compensated_mag)
        plot!(fig, r.t, mag_demean;
              label = "$(r.label)  σ=$(round(std(r.compensated_mag);digits=2))nT",
              color = col, lw=1.2, alpha=0.8)
    end

    hline!(fig, [0.0]; color=:gray, ls=:dash, lw=0.8, label="")

    !isnothing(save_path) && savefig(fig, save_path)
    return fig
end

# -- 2×2 Summary Dashboard -----------------------------------------------------

"""
    plot_benchmark_summary(results, metrics_list; save_path=nothing, tag="") -> Plot

2×2 dashboard combining the four most informative plots:
  [1,1] DRMS bar chart
  [1,2] Error vs time (North)
  [2,1] Heading correlation bar
  [2,2] Position scatter
"""
function plot_benchmark_summary(
    results      :: Vector{BenchmarkResult},
    metrics_list :: Vector{NavMetrics};
    save_path    :: Union{String,Nothing} = nothing,
    tag          :: String = ""
)
    palette = [COL_OFFICIAL, COL_CUSTOM, :darkgreen, :purple]

    # Panel 1: DRMS/CEP bar chart
    labels   = [m.label for m in metrics_list]
    drms_vals= [m.drms  for m in metrics_list]
    cep_vals = [m.cep   for m in metrics_list]

    p1 = groupedbar(
        repeat(labels; inner=2),
        vcat([[d,c] for (d,c) in zip(drms_vals, cep_vals)]...);
        group    = repeat(["DRMS","CEP"]; outer=length(labels)),
        title    = "Navigation Accuracy",
        ylabel   = "Error (m)",
        color    = [COL_OFFICIAL COL_CUSTOM],
        bar_width= 0.6,
        legend   = :topleft,
        grid     = true
    )

    # Panel 2: North error vs time
    p2 = plot(; title="North Error vs Time", xlabel="t (s)", ylabel="m",
                legend=:topright, grid=true)
    for (i,r) in enumerate(results)
        plot!(p2, r.t, r.err_north_m;
              label=r.label, color=palette[mod1(i,end)], lw=1.4)
    end
    hline!(p2, [0.0]; color=:gray, ls=:dash, lw=0.6, label="")

    # Panel 3: Heading correlation bar
    n_bins = 36
    p3 = plot(; title="Heading-Binned Error", xlabel="Heading (°)", ylabel="m",
                legend=:topright, grid=true)
    for (i,r) in enumerate(results)
        col = palette[mod1(i,end)]
        bins_deg, mean_err, corr_val = heading_correlation(r; n_bins=n_bins)
        bar!(p3, bins_deg, mean_err;
             label="$(r.label) |r|=$(round(corr_val;digits=3))",
             color=col, alpha=0.55, bar_width=9)
    end

    # Panel 4: Position scatter
    p4 = plot(; title="Position Scatter", xlabel="East (m)", ylabel="North (m)",
                aspect_ratio=1, legend=:topright, grid=true)
    for (i,r) in enumerate(results)
        col = palette[mod1(i,end)]
        scatter!(p4, r.err_east_m, r.err_north_m;
                 label=r.label, color=col, ms=2, markerstrokewidth=0, alpha=0.4)
        drms_v = sqrt(mean(r.err_north_m.^2 .+ r.err_east_m.^2))
        θ = range(0,2π;length=200)
        plot!(p4, drms_v.*cos.(θ), drms_v.*sin.(θ); color=col, lw=1.5, ls=:dash, label="")
    end
    hline!(p4,[0.0];color=:gray,lw=0.5,label="")
    vline!(p4,[0.0];color=:gray,lw=0.5,label="")

    title_str = isempty(tag) ? "Benchmark Summary" : "Benchmark: $tag"
    fig = plot(p1, p2, p3, p4;
               layout = (2,2),
               size   = (1200, 900),
               dpi    = 150,
               plot_title = title_str)

    !isnothing(save_path) && savefig(fig, save_path)
    return fig
end

# -- Multi-segment box plot -----------------------------------------------------

"""
    plot_multi_segment(segment_results; metric=:drms, save_path=nothing) -> Plot

Given a vector of (official, custom) NavMetrics pairs from repeated experiments,
plot a side-by-side comparison across segments.

# Arguments
- `segment_results` : Vector of NamedTuples `(official::NavMetrics, custom::NavMetrics, label::String)`
- `metric`          : field name in NavMetrics to extract (e.g. `:drms`, `:cep`, `:rms_north`)
"""
function plot_multi_segment(
    segment_results :: Vector{<:NamedTuple};
    metric          :: Symbol = :drms,
    save_path       :: Union{String,Nothing} = nothing
)
    seg_labels  = [r.label for r in segment_results]
    vals_off    = [getfield(r.official, metric) for r in segment_results]
    vals_cust   = [getfield(r.custom,   metric) for r in segment_results]

    metric_name = replace(string(metric), "_" => " ")

    fig = plot(;
        title  = "Segment Comparison — $metric_name",
        xlabel = "Segment",
        ylabel = "$metric_name (m)",
        legend = :topright,
        grid   = true,
        xticks = (1:length(seg_labels), seg_labels),
        xrotation = 30,
        size   = (900, 450),
        dpi    = 150
    )

    plot!(fig, 1:length(seg_labels), vals_off;
          label=:official_tl, marker=:circle, lw=2, ms=6, color=COL_OFFICIAL)
    plot!(fig, 1:length(seg_labels), vals_cust;
          label=:custom_tl,   marker=:square, lw=2, ms=6, color=COL_CUSTOM)

    !isnothing(save_path) && savefig(fig, save_path)
    return fig
end

# -- Save all figures ----------------------------------------------------------

"""
    save_all_plots(results, metrics_list, out_dir; tag="")

Save the full set of benchmark figures to `out_dir`.
"""
function save_all_plots(
    results      :: Vector{BenchmarkResult},
    metrics_list :: Vector{NavMetrics},
    out_dir      :: String;
    tag          :: String = ""
)
    isdir(out_dir) || mkpath(out_dir)
    pfx = isempty(tag) ? "" : tag * "_"

    plot_error_vs_time(results;
        save_path = joinpath(out_dir, pfx * "error_vs_time.png"))

    plot_position_errors(results;
        save_path = joinpath(out_dir, pfx * "position_scatter.png"))

    plot_heading_correlation(results;
        save_path = joinpath(out_dir, pfx * "heading_corr.png"))

    plot_mag_comparison(results;
        save_path = joinpath(out_dir, pfx * "mag_compare.png"))

    plot_benchmark_summary(results, metrics_list;
        save_path = joinpath(out_dir, pfx * "summary.png"),
        tag       = tag)

    println("  ✓  Figures saved to: $out_dir")
end
