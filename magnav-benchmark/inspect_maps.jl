#!/usr/bin/env julia
"""
    inspect_maps.jl
    
Inspect the MapS struct from MagNav package.
"""

using MagNav

println("\n" * "="^60)
println("  Inspecting MapS struct from MagNav")
println("="^60 * "\n")

# Try to get the docstring for MapS
try
    println("MapS Docstring:")
    println("-" * 60)
    @doc MapS
    println()
catch e
    println("Could not retrieve MapS docstring: $e\n")
end

# Try to print the struct type info
try
    println("MapS Type Info:")
    println("-" * 60)
    println(MapS)
    println()
catch e
    println("Could not inspect MapS type: $e\n")
end

# Try to find constructors
try
    println("MapS Methods (constructors):")
    println("-" * 60)
    methods(MapS)
    println()
catch e
    println("Could not list methods: $e\n")
end

# Check what modules export MapS
try
    println("Exported from modules:")
    println("-" * 60)
    # This is indirect inspection
    println("MapS is available from: MagNav")
    println()
catch e
    println("Could not determine exporting modules: $e\n")
end

println("="^60)
println("Inspection complete")
println("="^60 * "\n")
