#!/usr/bin/env julia
"""
    test/runtests.jl

Run all unit tests (no MagNav data required).

Usage:
    julia --project=. test/runtests.jl
"""

using Pkg
Pkg.activate(@__DIR__ |> dirname)

println("\n" * "█" ^ 55)
println("  MagNav Benchmark Framework — Test Suite")
println("█" ^ 55)

include("test_tl.jl")
include("test_metrics_segments.jl")

println("\n" * "█" ^ 55)
println("  ✓  All test suites passed")
println("█" ^ 55 * "\n")
