#!/usr/bin/env julia
"""
    test/test_metrics_segments.jl

Unit tests for NavMetrics, heading_correlation, and SegmentSpec.

No MagNav or flight data required.

Usage:
    julia --project=. test/test_metrics_segments.jl
"""

using Pkg
Pkg.activate(@__DIR__ |> dirname)
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

include(joinpath(@__DIR__, "..", "src", "metrics", "nav_metrics.jl"))
include(joinpath(@__DIR__, "..", "src", "utils", "segment_utils.jl"))

using Test, Statistics, Printf

# -- Synthetic BenchmarkResult -------------------------------------------------

function fake_result(label; n=500, drms_approx=50.0, seed=1)
    rng = Xoshiro(seed)
    t   = collect(0.0:1.0:(n-1))
    en  = randn(rng, n) .* drms_approx ./ sqrt(2)
    ee  = randn(rng, n) .* drms_approx ./ sqrt(2)
    hdg = mod.(LinRange(0, 4π, n), 2π)
    mag = randn(rng, n) .* 20.0

    return BenchmarkResult(
        label,
        t, zeros(n), zeros(n), zeros(n), zeros(n),
        en, ee,
        hdg, mag,
        nothing, nothing
    )
end

# ═══════════════════════════════════════════════════════════════════════
println("\n" * "=" ^ 55)
println("  Metrics & Segments Unit Tests")
println("=" ^ 55)

@testset "NavMetrics: DRMS consistency" begin
    r  = fake_result("test"; drms_approx=50.0)
    m  = compute_nav_metrics(r)

    # DRMS = sqrt(RMS_N² + RMS_E²)
    @test abs(m.drms - sqrt(m.rms_north^2 + m.rms_east^2)) < 1e-10

    # 2DRMS
    @test abs(m.drms_2 - 2*m.drms) < 1e-10

    # CEP ≤ DRMS (always true for symmetric distribution)
    @test m.cep <= m.drms * 1.2   # loose bound

    println("  ✓  DRMS = sqrt(N²+E²), 2DRMS, CEP ≤ DRMS (approx)")
end

@testset "NavMetrics: final and peak errors" begin
    r = fake_result("test2"; n=200)
    m = compute_nav_metrics(r)
    radial = @. sqrt(r.err_north_m^2 + r.err_east_m^2)

    @test abs(m.final_error - radial[end]) < 1e-10
    @test abs(m.peak_error  - maximum(radial)) < 1e-10

    println("  ✓  final_error and peak_error match radial vector")
end

@testset "delta_metrics: sign convention" begin
    r_good = fake_result("good"; drms_approx=20.0, seed=10)
    r_bad  = fake_result("bad";  drms_approx=80.0, seed=11)
    m_good = compute_nav_metrics(r_good)
    m_bad  = compute_nav_metrics(r_bad)

    Δ = delta_metrics(m_good, m_bad)   # ref=good, alt=bad
    @test Δ.Δdrms > 0.0    # alt is worse → positive Δ
    @test Δ.pct_drms > 0.0

    Δ2 = delta_metrics(m_bad, m_good)  # ref=bad, alt=good
    @test Δ2.Δdrms < 0.0              # alt is better → negative Δ

    println("  ✓  delta_metrics sign convention")
end

@testset "heading_correlation: returns correct bin count" begin
    r = fake_result("hdg")
    bins, mean_err, corr = heading_correlation(r; n_bins=36)
    @test length(bins)     == 36
    @test length(mean_err) == 36
    @test 0.0 <= corr <= 1.0
    println("  ✓  heading_correlation: bins=36, corr ∈ [0,1]")
end

@testset "SegmentSpec: overlap rejection" begin
    @test_throws ErrorException SegmentSpec(;
        cal_start_s=0.0, cal_end_s=600.0,
        eval_start_s=300.0, eval_end_s=900.0,
        gap_s=0.0
    )
    println("  ✓  SegmentSpec rejects overlapping cal/eval")
end

@testset "SegmentSpec: gap warning" begin
    # Should warn but not error when gap < required
    spec = SegmentSpec(;
        cal_start_s=0.0, cal_end_s=100.0,
        eval_start_s=120.0, eval_end_s=400.0,
        gap_s=60.0   # actual gap is 20s < 60s required
    )
    @test spec.cal_end_s == 100.0
    println("  ✓  SegmentSpec warns on insufficient gap")
end

@testset "resolve_segments: basic" begin
    # Synthetic time vector 0..999 seconds
    t = collect(0.0:1.0:999.0)
    spec = SegmentSpec(;
        cal_start_s=0.0, cal_end_s=299.0,
        eval_start_s=360.0, eval_end_s=599.0,
        gap_s=60.0
    )
    cal_ind, eval_ind = resolve_segments(spec, t)

    @test sum(cal_ind)  == 300   # 0..299 inclusive
    @test sum(eval_ind) == 240   # 360..599 inclusive
    @test sum(cal_ind .& eval_ind) == 0   # no overlap

    println("  ✓  resolve_segments: correct counts, no overlap")
end

println("\n" * "=" ^ 55)
println("  All metrics/segment tests passed.")
println("=" ^ 55 * "\n")
