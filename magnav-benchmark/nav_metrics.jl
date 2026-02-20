"""
    nav_metrics.jl

Navigation-level metrics for benchmarking TL compensation quality.

All metrics operate on `BenchmarkResult` structs and return typed `NavMetrics`.

Metrics provided
----------------
| Name                    | Symbol      | Unit    |
|-------------------------|-------------|---------|
| RMS North error         | rms_north   | m       |
| RMS East error          | rms_east    | m       |
| Distance RMS (DRMS)     | drms        | m       |
| 2DRMS (95th percentile) | drms_2      | m       |
| CEP (50th pct circle)   | cep         | m       |
| Final position error    | final_error | m       |
| Peak error              | peak_error  | m       |
| Mean heading correlation| hdg_corr    | -       |
| Residual RMS (mag)      | mag_rms     | nT      |
"""

# -- Result type ---------------------------------------------------------------

"""
    NavMetrics

Typed container for all navigation benchmark metrics.
"""
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

# -- Core computation ----------------------------------------------------------

"""
    compute_nav_metrics(result::BenchmarkResult) -> NavMetrics

Compute all navigation metrics for one benchmark run.
"""
function compute_nav_metrics(result::BenchmarkResult)
    en = result.err_north_m
    ee = result.err_east_m

    radial    = @. sqrt(en^2 + ee^2)

    rms_n     = sqrt(mean(en.^2))
    rms_e     = sqrt(mean(ee.^2))
    drms_val  = sqrt(rms_n^2 + rms_e^2)
    drms2_val = 2 * drms_val
    cep_val   = quantile(radial, 0.50)
    final_val = radial[end]
    peak_val  = maximum(radial)

    # Heading correlation: correlation between compensated mag and heading
    # A well-compensated signal should have near-zero correlation
    hdg_corr_val = abs(cor(result.compensated_mag, result.heading_rad))

    # Residual RMS of the compensated scalar (after mean removal)
    mag_rms_val  = sqrt(mean((result.compensated_mag .- mean(result.compensated_mag)).^2))

    return NavMetrics(
        result.label,
        rms_n, rms_e, drms_val, drms2_val,
        cep_val, final_val, peak_val,
        hdg_corr_val, mag_rms_val
    )
end

"""
    compute_nav_metrics(results::Vector{BenchmarkResult}) -> Vector{NavMetrics}

Batch version: compute metrics for a list of runs.
"""
compute_nav_metrics(results::Vector{<:Any}) = compute_nav_metrics.(results)

# -- Heading correlation -------------------------------------------------------

"""
    heading_correlation(result::BenchmarkResult; n_bins=36) -> (bins, mean_err, corr)

Compute heading-binned mean position error and linear correlation.

Returns:
- `bins`     : heading bin centres (deg)
- `mean_err` : mean radial error per bin (m)
- `corr`     : Pearson |r| between heading (deg) and radial error
"""
function heading_correlation(result::BenchmarkResult; n_bins::Int=36)
    hdg_deg  = rad2deg.(result.heading_rad)
    radial   = @. sqrt(result.err_north_m^2 + result.err_east_m^2)

    edges    = range(0.0, 360.0; length=n_bins+1)
    centres  = collect(edges[1:end-1] .+ step(edges)/2)
    mean_err = zeros(n_bins)

    for (i, (lo, hi)) in enumerate(zip(edges[1:end-1], edges[2:end]))
        mask = (hdg_deg .>= lo) .& (hdg_deg .< hi)
        mean_err[i] = sum(mask) > 0 ? mean(radial[mask]) : NaN
    end

    corr_val = abs(cor(hdg_deg, radial))

    return (bins=centres, mean_err=mean_err, corr=corr_val)
end

# -- Comparison table ---------------------------------------------------------

"""
    metrics_table(metrics_list::Vector{NavMetrics}) -> DataFrame

Convert a list of NavMetrics to a formatted comparison DataFrame.
Requires DataFrames.jl (imported in module).
"""
function metrics_table(metrics_list::Vector{NavMetrics})
    df = DataFrame(
        Method      = [m.label      for m in metrics_list],
        RMS_N_m     = [m.rms_north  for m in metrics_list],
        RMS_E_m     = [m.rms_east   for m in metrics_list],
        DRMS_m      = [m.drms       for m in metrics_list],
        DRMS2_m     = [m.drms_2     for m in metrics_list],
        CEP_m       = [m.cep        for m in metrics_list],
        Final_m     = [m.final_error for m in metrics_list],
        Peak_m      = [m.peak_error  for m in metrics_list],
        Hdg_Corr    = [m.hdg_corr   for m in metrics_list],
        MagRMS_nT   = [m.mag_rms    for m in metrics_list],
    )
    return df
end

"""
    print_metrics_table(metrics_list)

Print a formatted ASCII comparison table to stdout.
"""
function print_metrics_table(metrics_list::Vector{NavMetrics})
    hdr = @sprintf("%-20s %8s %8s %8s %8s %8s %8s %8s %8s %8s", "Method", "RMS_N", "RMS_E", "DRMS", "2DRMS", "CEP", "Final", "Peak", "|r_hdg|", "MagRMS")
    sep = "-" ^ length(hdr)

    println(sep)
    println(hdr)
    println(sep)
    for m in metrics_list
        @printf("%-20s %7.2fm %7.2fm %7.2fm %7.2fm %7.2fm %7.2fm %7.2fm %8.4f %7.2fnT\n", m.label, m.rms_north, m.rms_east, m.drms, m.drms_2, m.cep, m.final_error, m.peak_error, m.hdg_corr, m.mag_rms)
    end
    println(sep)
end

"""
    delta_metrics(ref::NavMetrics, alt::NavMetrics) -> NamedTuple

Compute signed difference (alt - ref) for all scalar metrics.
Positive = alt is worse; negative = alt is better.
"""
function delta_metrics(ref::NavMetrics, alt::NavMetrics)
    return (
        Δrms_north   = alt.rms_north   - ref.rms_north,
        Δrms_east    = alt.rms_east    - ref.rms_east,
        Δdrms        = alt.drms        - ref.drms,
        Δcep         = alt.cep         - ref.cep,
        Δfinal_error = alt.final_error - ref.final_error,
        Δpeak_error  = alt.peak_error  - ref.peak_error,
        Δhdg_corr    = alt.hdg_corr    - ref.hdg_corr,
        Δmag_rms     = alt.mag_rms     - ref.mag_rms,
        pct_drms     = 100 * (alt.drms - ref.drms) / ref.drms,
    )
end
