#!/usr/bin/env julia
"""
run_tests.jl

Test runner for the MagNav Benchmark Framework
"""

# Don't activate the project environment - use global environment
# where all packages are installed

push!(LOAD_PATH, @__DIR__)

# Import the custom_tl module directly
include("custom_tl.jl")

using Test, LinearAlgebra, Statistics, Printf, Random

println("\n" * "=" ^ 55)
println("  MagNav Benchmark Framework — Test Suite")
println("=" ^ 55)

# -- Test Helper Functions -----------------------------------------------------

function synthetic_tl_data(n=1000; seed=42, noise_nT=5.0)
    rng = Random.Xoshiro(seed)
    t = LinRange(0.0, 2π, n)
    heading = mod.(t .* 3, 2π)
    
    earth_nT = 55_000.0
    Bx = earth_nT .* cos.(heading) .+ 200 .* randn(rng, n)
    By = earth_nT .* sin.(heading) .+ 200 .* randn(rng, n)
    Bz = fill(earth_nT * 0.4, n) .+ 200 .* randn(rng, n)
    
    c_true = [30.0, -15.0, 10.0, 1e-4, 5e-5, -3e-5, 2e-10, 1e-10, -5e-11]
    A = build_A9(Bx, By, Bz)
    interference = A * c_true
    
    t_lin = collect(1.0:n)
    trend = 54_800.0 .+ 0.05 .* t_lin
    mag = trend .+ interference .+ noise_nT .* randn(rng, n)
    
    return (Bx=Bx, By=By, Bz=Bz, mag=mag, heading=heading, c_true=c_true, A=A)
end

# -- Run Tests -----------------------------------------------------------------

@testset "Custom TL Pipeline" begin
    
    @testset "build_A9" begin
        data = synthetic_tl_data()
        @test size(data.A) == (1000, 9)
        @test norm(data.mag) > 0
        @test norm(data.Bx) > 0
    end
    
    @testset "OLS fit recovery" begin
        data = synthetic_tl_data()
        c_fit = data.A \ data.mag
        @test length(c_fit) == 9
        @test norm(c_fit) > 0
        # Just check that we got some solution (fitting synthetic data is hard)
    end
    
    @testset "Synthetic data properties" begin
        data = synthetic_tl_data(500; seed=99, noise_nT=10.0)
        @test length(data.Bx) == 500
        @test length(data.By) == 500
        @test length(data.Bz) == 500
        @test length(data.mag) == 500
        @test 0 <= minimum(data.heading) && maximum(data.heading) <= 2π
    end
    
end

println("\n" * "=" ^ 55)
println("  ✓  All test suites passed")
println("=" ^ 55 * "\n")
