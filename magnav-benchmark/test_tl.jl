#!/usr/bin/env julia
"""
    test/test_tl.jl

Unit tests for the custom TL pipeline.

Tests are synthetic (no MagNav data required) and run in < 2 seconds.

Usage:
    julia --project=. test/test_tl.jl
"""

using Pkg
Pkg.activate(@__DIR__ |> dirname)
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

# Use only the TL modules (no MagNav needed for unit tests)
include(joinpath(@__DIR__, "..", "src", "tl", "custom_tl.jl"))

using Test, LinearAlgebra, Statistics, Printf

# -- Helpers -------------------------------------------------------------------

function synthetic_tl_data(n=1000; seed=42, noise_nT=5.0)
    rng = Xoshiro(seed)

    # Simulate aircraft flying a cloverleaf: heading cycles through 360°
    t       = LinRange(0.0, 2π, n)
    heading = mod.(t .* 3, 2π)     # 3 full rotations

    # Synthetic fluxgate (body frame)
    earth_nT = 55_000.0
    Bx = earth_nT .* cos.(heading) .+ 200 .* randn(rng, n)
    By = earth_nT .* sin.(heading) .+ 200 .* randn(rng, n)
    Bz = fill(earth_nT * 0.4, n)  .+ 200 .* randn(rng, n)

    # True TL coefficients (9-term)
    c_true = [30.0, -15.0, 10.0,   # permanent
              1e-4,  5e-5, -3e-5,   # induced linear
              2e-10, 1e-10, -5e-11] # induced quadratic

    A = build_A9(Bx, By, Bz)
    interference = A * c_true

    # Noisy scalar measurement
    t_lin  = collect(1.0:n)
    trend  = 54_800.0 .+ 0.05 .* t_lin
    mag    = trend .+ interference .+ noise_nT .* randn(rng, n)

    return (Bx=Bx, By=By, Bz=Bz, mag=mag, heading=heading, c_true=c_true, A=A)
end

# ═══════════════════════════════════════════════════════════════════════
println("\n" * "=" ^ 55)
println("  Custom TL Unit Tests")
println("=" ^ 55)

@testset "build_A9" begin
    n = 100
    Bx = randn(n); By = randn(n); Bz = randn(n)
    A = build_A9(Bx, By, Bz)

    @test size(A) == (n, 9)
    @test !any(isnan, A)
    @test !any(isinf, A)

    # Unit direction cosines: cols 1-3 should have |row| ≤ 1
    norms_perm = @. sqrt(A[:,1]^2 + A[:,2]^2 + A[:,3]^2)
    @test all(abs.(norms_perm .- 1.0) .< 1e-10)

    println("  ✓  build_A9: shape, no NaN, unit direction cosines")
end

@testset "fit_custom_tl: recovers synthetic coefficients" begin
    d = synthetic_tl_data(2000; noise_nT=1.0)

    model = CustomTLModel()
    fit_custom_tl!(model, d.Bx, d.By, d.Bz, d.mag, d.heading;
                   lambda=0.0, detrend=true, verbose=false)

    @test model.fitted
    @test length(model.coef) == 9

    # Coefficient recovery: within 20% of each true coef (noise makes this fuzzy)
    for (i, (ĉ, c)) in enumerate(zip(model.coef, d.c_true))
        err_pct = abs(ĉ - c) / (abs(c) + 1e-12) * 100
        @test err_pct < 20.0  "Coef $i: got $ĉ, expected $c (err=$(round(err_pct;digits=1))%)"
    end

    # Variance reduction should be positive
    @test model.var_reduction > 0.0

    # Heading correlation should drop
    @test model.heading_corr_post < model.heading_corr_pre

    println("  ✓  fit_custom_tl: coefficients, VR, heading correlation")
end

@testset "fit_custom_tl: variance reduction increases with lower noise" begin
    vr_noisy = let d = synthetic_tl_data(2000; noise_nT=50.0)
        m = CustomTLModel()
        fit_custom_tl!(m, d.Bx, d.By, d.Bz, d.mag, d.heading;
                       detrend=true, verbose=false)
        m.var_reduction
    end

    vr_clean = let d = synthetic_tl_data(2000; noise_nT=1.0)
        m = CustomTLModel()
        fit_custom_tl!(m, d.Bx, d.By, d.Bz, d.mag, d.heading;
                       detrend=true, verbose=false)
        m.var_reduction
    end

    @test vr_clean > vr_noisy
    println("  ✓  VR monotone in noise level: vr_clean=$(@sprintf "%.3f" vr_clean) > vr_noisy=$(@sprintf "%.3f" vr_noisy)")
end

@testset "compensate_custom_tl: residuals reduce interference" begin
    d     = synthetic_tl_data(2000; noise_nT=2.0)
    model = CustomTLModel()
    fit_custom_tl!(model, d.Bx, d.By, d.Bz, d.mag, d.heading;
                   detrend=true, verbose=false)

    comp = compensate_custom_tl(model, d.Bx, d.By, d.Bz, d.mag; detrend=true)

    @test length(comp) == length(d.mag)
    @test !any(isnan, comp)

    # RMS of compensated should be < RMS of raw (after detrend)
    t  = collect(1.0:length(d.mag))
    tc = [ones(length(t)) t] \ d.mag
    y_dt = d.mag .- [ones(length(t)) t] * tc
    @test sqrt(mean(comp.^2)) < sqrt(mean(y_dt.^2))

    println("  ✓  compensate_custom_tl: reduces signal RMS")
end

@testset "ridge regression: lambda > 0 reduces coefficient norm" begin
    d = synthetic_tl_data(500; noise_nT=10.0)
    m0 = CustomTLModel()
    m1 = CustomTLModel()
    fit_custom_tl!(m0, d.Bx, d.By, d.Bz, d.mag, d.heading;
                   lambda=0.0, detrend=true, verbose=false)
    fit_custom_tl!(m1, d.Bx, d.By, d.Bz, d.mag, d.heading;
                   lambda=1e3, detrend=true, verbose=false)

    @test norm(m1.coef) < norm(m0.coef)
    println("  ✓  ridge: λ>0 shrinks coefficient norm")
end

@testset "CustomTLModel: unfitted guard" begin
    m = CustomTLModel()
    @test_throws AssertionError compensate_custom_tl(m, [1.0], [1.0], [1.0], [1.0])
    println("  ✓  unfitted guard raises AssertionError")
end

println("\n" * "=" ^ 55)
println("  All tests passed.")
println("=" ^ 55 * "\n")
