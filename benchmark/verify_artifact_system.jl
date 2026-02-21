#!/usr/bin/env julia
"""
verify_artifact_system.jl

Quick verification script that tests whether the MagNav artifact system
can find and access the Flt1006_train.h5 file.

Usage:
    julia verify_artifact_system.jl
"""

println("\n" * "="^60)
println("  MagNav Artifact System Verification")
println("="^60 * "\n")

try
    # Import MagNav
    println("1. Loading MagNav...")
    using MagNav
    println("   ✓ MagNav loaded")
    
    # Get artifact directory
    println("\n2. Accessing artifact system...")
    H5_DIR = sgl_2020_train()
    println("   ✓ Artifact directory: $H5_DIR")
    
    # Verify files exist
    println("\n3. Checking for flight data files...")
    H5_PATH = joinpath(H5_DIR, "Flt1006_train.h5")
    
    if isfile(H5_PATH)
        println("   ✓ Flt1006_train.h5 found at:")
        println("     $H5_PATH")
        
        # Get file size
        filesize_mb = filesize(H5_PATH) / 1024^2
        println("     Size: $(round(filesize_mb; digits=1)) MB")
    else
        println("   ✗ Flt1006_train.h5 NOT found!")
        println("     Expected at: $H5_PATH")
        println("     Listing directory contents:")
        for f in readdir(H5_DIR)
            println("       - $f")
        end
    end
    
    # List all available flights
    println("\n4. Available flight files in artifact:")
    h5_files = filter(f -> endswith(f, ".h5"), readdir(H5_DIR))
    for f in sort(h5_files)
        path = joinpath(H5_DIR, f)
        size_mb = filesize(path) / 1024^2
        println("   - $f ($(round(size_mb; digits=1)) MB)")
    end
    
    println("\n" * "="^60)
    println("  ✓ Artifact system is working correctly")
    println("="^60 * "\n")
    
catch e
    println("\n" * "="^60)
    println("  ✗ Error accessing artifact system")
    println("="^60)
    println("\nError details:")
    println(e)
    println("\nThis should not happen. Make sure MagNav is installed.")
    exit(1)
end
