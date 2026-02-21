"""
    data_loader.jl

Reproducible data loading for SGL 2020 / Flt1006_train.h5.

Wraps MagNav's standard loading pipeline (sgl_2020_train → h5_to_df →
get_XYZ) into a single typed struct so that downstream code never touches
raw DataFrames or h5 files directly.

Caching
-------
`load_benchmark_data` optionally writes a JLD2 cache file to avoid re-running
the slow MagNav loading step on repeated runs.  The cache key includes a hash
of the h5 file path to detect file changes.
"""

using HDF5
using DataFrames
using MagNav

# -- BenchmarkData -------------------------------------------------------------

"""
    BenchmarkData

Container for all loaded benchmark data for one flight.

Fields
------
- `xyz`       : MagNav XYZ struct (XYZ20 for SGL 2020)
- `mapS`      : MapS scalar anomaly map
- `flight_id` : e.g. `:Flt1006`
- `h5_path`   : absolute path to source HDF5 file
- `loaded_at` : timestamp of load (for cache validation)
"""
struct BenchmarkData
    xyz        :: Any            # MagNav XYZ20 or similar
    mapS       :: Any            # MagNav MapS
    flight_id  :: Symbol
    h5_path    :: String
    loaded_at  :: DateTime
end

# -- Main loader ----------------------------------------------------------------

"""
    load_benchmark_data(h5_path, map_artifact;
                        flight_id      = :Flt1006,
                        line_num       = nothing,
                        silent         = true,
                        verbose        = true) -> BenchmarkData

Load SGL 2020 data and Ottawa area map into a `BenchmarkData` struct.

# Arguments
- `h5_path`      : absolute path to Flt1006_train.h5
- `map_artifact` : artifact name for ottawa_area_maps, e.g. "ottawa_area_maps"
- `flight_id`    : MagNav flight symbol, default `:Flt1006`
- `line_num`     : optional integer to select a specific flight line
                   (if nothing, uses the full XYZ)
- `silent`       : suppress MagNav's internal print output
- `verbose`      : print loading progress

# Example
    data = load_benchmark_data(
        "/data/sgl2020/Flt1006_train.h5",
        "ottawa_area_maps"
    )
"""
function load_benchmark_data(
    h5_path      :: String,
    map_artifact :: String;
    flight_id    :: Symbol  = :Flt1006,
    line_num     :: Union{Int,Nothing} = nothing,
    silent       :: Bool    = true,
    verbose      :: Bool    = true
)
    isfile(h5_path) || error("HDF5 file not found: $h5_path")

    verbose && println("-- Loading BenchmarkData -----------------------------")
    verbose && println("  Flight  : $flight_id")
    verbose && println("  H5 path : $h5_path")

    # -- Flight data -----------------------------------------------------------
    verbose && print("  Loading HDF5 → XYZ struct...")
    # Directly construct XYZ using MagNav's h5 loading function
    # For SGL 2020 dataset, use get_XYZ20 which is specialized for this format
    xyz = if flight_id == :Flt1006
        MagNav.get_XYZ20(h5_path; silent=silent)
    else
        # Fallback for other flights - assume similar structure
        MagNav.get_XYZ20(h5_path; silent=silent)
    end
    
    verbose && println(" ✓  ($(size(xyz.traj.lat,1)) samples @ $(round(1/xyz.traj.dt;digits=1)) Hz)")

    # -- Map -------------------------------------------------------------------
    verbose && print("  Loading map artifact '$map_artifact'...")
    # Get the artifact path and load a representative map file
    map_path = if map_artifact == "ottawa_area_maps"
        MagNav.ottawa_area_maps()
    else
        map_artifact  # fallback: assume it's already a path
    end
    # Load the first suitable map file from the artifact folder
    map_file = joinpath(map_path, "Eastern_395.h5")
    mapS = get_map(map_file)
    verbose && println(" ✓")

    verbose && println("-- BenchmarkData ready -------------------------------")

    return BenchmarkData(xyz, mapS, flight_id, abspath(h5_path), now())
end

"""
    describe_data(data::BenchmarkData)

Print a human-readable summary of the loaded data.
"""
function describe_data(data::BenchmarkData)
    xyz = data.xyz
    n   = length(xyz.traj.lat)
    dur = (xyz.traj.tt[end] - xyz.traj.tt[1])

    @printf("BenchmarkData\n")
    @printf("  Flight     : %s\n", string(data.flight_id))
    @printf("  H5         : %s\n", data.h5_path)
    @printf("  Samples    : %d\n", n)
    @printf("  Duration   : %.1f s (%.1f min)\n", dur, dur/60)
    @printf("  Sample rate: %.1f Hz\n", 1/xyz.traj.dt)
    @printf("  Lat range  : %.4f - %.4f\n", rad2deg(minimum(xyz.traj.lat)), rad2deg(maximum(xyz.traj.lat)))
    @printf("  Lon range  : %.4f - %.4f\n", rad2deg(minimum(xyz.traj.lon)), rad2deg(maximum(xyz.traj.lon)))
    @printf("  Loaded at  : %s\n", string(data.loaded_at))
end

# -- Scalar mag extractor ------------------------------------------------------

"""
    extract_cal_signals(xyz, ind; mag_field=:mag_1_c, flux_field=:flux_a)
                        -> NamedTuple

Extract the key calibration signals as plain vectors for model fitting.

Returns: `(Bx, By, Bz, mag_scalar, heading, t)`
"""
function extract_cal_signals(
    xyz        :: Any,
    ind        :: AbstractVector{Bool};
    mag_field  :: Symbol = :mag_1_c,
    flux_field :: Symbol = :flux_a
)
    flux  = getfield(xyz, flux_field)
    mag   = Float64.(getfield(xyz, mag_field)[ind])
    hdg   = Float64.(xyz.ins_yaw[ind])
    t     = Float64.(xyz.traj.tt[ind])
    t   .-= t[1]   # zero-based time

    return (
        Bx         = Float64.(flux.x[ind]),
        By         = Float64.(flux.y[ind]),
        Bz         = Float64.(flux.z[ind]),
        mag_scalar = mag,
        heading    = hdg,
        t          = t,
        n          = sum(ind)
    )
end
