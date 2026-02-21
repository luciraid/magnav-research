#!/usr/bin/env julia
"""
export_ottawa_map.jl — Export Ottawa area magnetic anomaly map.

Run once before launching the C++ plugin or Python bridge:
    julia gazebo/plugins/export_ottawa_map.jl

Outputs (in gazebo/plugins/data/):
  ottawa_map.npz   — NumPy format for the Python magnetic_map_bridge.py
  ottawa_map.bin   — Raw binary format for the C++ MagneticMapPlugin

Binary format (ottawa_map.bin) — little-endian:
    uint32_t  nlat          number of latitude  grid points
    uint32_t  nlon          number of longitude grid points
    float64_t lat0_rad      latitude  of first grid point  (radians)
    float64_t dlat_rad      latitude  grid step             (radians)
    float64_t lon0_rad      longitude of first grid point  (radians)
    float64_t dlon_rad      longitude grid step             (radians)
    float64_t alt_m         survey altitude above WGS-84   (metres)
    float32_t data[nlat*nlon]  row-major anomaly values    (nT, NaN=masked)

Map: Eastern_395.h5 — Eastern Ontario survey flown by SGL, upward-continued
     to 395 m above WGS-84. Values are total-field magnetic anomaly in nT.
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", "MagNav.jl"), io=devnull)
using MagNav

# ── Output path ───────────────────────────────────────────────────────────────
out_dir  = joinpath(@__DIR__, "data")
out_path = joinpath(out_dir, "ottawa_map.npz")
mkpath(out_dir)

# ── Load map ──────────────────────────────────────────────────────────────────
map_dir  = ottawa_area_maps()
map_file = joinpath(map_dir, "Eastern_395.h5")
@info "Loading map from: $map_file"
m = get_map(map_file)

# Validate
@assert eltype(m.yy) == Float64 "Expected Float64 lat vector"
@assert eltype(m.xx) == Float64 "Expected Float64 lon vector"
@assert size(m.map) == (length(m.yy), length(m.xx)) "Unexpected map dimensions"

# ── Apply mask: fill outside-survey cells with NaN ───────────────────────────
map_nT = Float64.(m.map)
map_nT[.!m.mask] .= NaN

# ── Write NPZ via npzwrite ────────────────────────────────────────────────────
# Use NPZ.jl (add if needed: Pkg.add("NPZ"))
try
    using NPZ
catch
    @info "Adding NPZ.jl..."
    Pkg.add("NPZ")
    using NPZ
end

npzwrite(out_path, Dict(
    "lat_rad" => Float64.(m.yy),
    "lon_rad" => Float64.(m.xx),
    "map_nT"  => Float64.(map_nT),
    "alt_m"   => Float64[m.alt],
))

# ── Write binary format for C++ plugin ───────────────────────────────────────
# Header: nlat(u32), nlon(u32), lat0(f64), dlat(f64), lon0(f64), dlon(f64), alt(f64)
# Data:   float32[nlat*nlon] row-major  (NaN for masked cells)
bin_path = joinpath(out_dir, "ottawa_map.bin")
let
    nlat   = UInt32(length(m.yy))
    nlon   = UInt32(length(m.xx))
    lat0   = Float64(m.yy[1])
    dlat   = Float64(m.yy[2] - m.yy[1])
    lon0   = Float64(m.xx[1])
    dlon   = Float64(m.xx[2] - m.xx[1])
    alt    = Float64(m.alt)
    data32 = Float32.(map_nT)   # NaN propagates through the cast

    open(bin_path, "w") do io
        write(io, nlat); write(io, nlon)
        write(io, lat0); write(io, dlat)
        write(io, lon0); write(io, dlon)
        write(io, alt)
        write(io, data32)
    end
    @info "Binary map written: $bin_path ($(round(filesize(bin_path)/1e6, digits=1)) MB)"
end

# ── Summary ───────────────────────────────────────────────────────────────────
lat_deg = extrema(rad2deg.(m.yy))
lon_deg = extrema(rad2deg.(m.xx))
println("=" ^ 60)
println("  Ottawa map exported successfully")
println("  NPZ    : $out_path")
println("  Binary : $bin_path")
println("  Shape  : $(size(map_nT)) (lat × lon)")
println("  Lat    : $(round.(lat_deg, digits=4)) deg")
println("  Lon    : $(round.(lon_deg, digits=4)) deg")
println("  Alt    : $(m.alt) m")
val_range = extrema(filter(!isnan, map_nT))
println("  Values : $(round.(val_range, digits=1)) nT")
println("  Valid  : $(sum(m.mask)) / $(length(m.mask)) cells")
println("=" ^ 60)
println("Python bridge:  python3 gazebo/plugins/magnetic_map_bridge.py")
println("C++ plugin:     build with  cmake gazebo/plugins && make")
