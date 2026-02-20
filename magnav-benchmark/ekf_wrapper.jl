"""
    ekf_wrapper.jl

Thin, side-effect-free wrapper around MagNav's `run_filt` for benchmarking.

Runs the EKF with a given compensated-scalar injection, captures the full
position output, and returns a typed `BenchmarkResult` for downstream metrics
and plotting.

Design principle: this file contains NO compensation logic.  It only:
  1. Injects the pre-computed signal into a copy of XYZ
  2. Calls `run_filt`
  3. Wraps the output

All compensation decisions (which model, which segments) are made by the
caller (the experiment script).
"""

# -- Result type ---------------------------------------------------------------

"""
    BenchmarkResult

Typed container for a single EKF navigation run.

Fields
------
- `label`         : human-readable run name, e.g. "official_tl" or "custom_tl"
- `t`             : time vector (seconds from epoch, length N)
- `lat_est`       : estimated latitude (rad), N
- `lon_est`       : estimated longitude (rad), N
- `lat_true`      : truth latitude (rad), N
- `lon_true`      : truth longitude (rad), N
- `err_north_m`   : northing error in metres, N
- `err_east_m`    : easting error in metres, N
- `heading_rad`   : true heading (rad), N — for correlation diagnostics
- `compensated_mag`: the scalar signal fed into the EKF (nT), N
- `filt_out`      : raw MagNav FILTout object (for advanced inspection)
- `config`        : ExperimentConfig snapshot (for reproducibility)
"""
struct BenchmarkResult
    label           :: String
    t               :: Vector{Float64}
    lat_est         :: Vector{Float64}
    lon_est         :: Vector{Float64}
    lat_true        :: Vector{Float64}
    lon_true        :: Vector{Float64}
    err_north_m     :: Vector{Float64}
    err_east_m      :: Vector{Float64}
    heading_rad     :: Vector{Float64}
    compensated_mag :: Vector{Float64}
    filt_out        :: Any          # MagNav FILTout (don't constrain type)
    config          :: Any          # ExperimentConfig
end

# -- Main entry point ----------------------------------------------------------

"""
    run_ekf_benchmark(label, xyz_injected, xyz_orig, mapS, ind;
                      config=default_config(), verbose=true) -> BenchmarkResult

Run MagNav EKF using `xyz_injected` (whose mag field has been replaced) and
return a `BenchmarkResult`.

# Arguments
- `label`         : name for this run
- `xyz_injected`  : copy of XYZ with the benchmark compensation applied
- `xyz_orig`      : original XYZ (used to extract truth/heading; not modified)
- `mapS`          : MagNav MapS struct (the Ottawa area scalar map)
- `ind`           : BitVector selecting the evaluation segment
- `config`        : ExperimentConfig (filter noise parameters, etc.)
- `verbose`       : print progress
"""
function run_ekf_benchmark(
    label        :: String,
    xyz_injected :: Any,
    xyz_orig     :: Any,
    mapS         :: Any,
    ind          :: AbstractVector{Bool};
    config       :: Any = nothing,
    verbose      :: Bool = true
)
    cfg = isnothing(config) ? default_config() : config

    verbose && @printf "\n▶ Running EKF: %s (%d samples)\n" label sum(ind)

    # -- Build map interpolant -------------------------------------------------
    map_vals, itp_mapS = get_map_val(mapS,
                                     xyz_injected.traj.lat[ind],
                                     xyz_injected.traj.lon[ind],
                                     xyz_injected.traj.alt[ind];
                                     return_itp=true)

    # -- Run EKF ---------------------------------------------------------------
    filt_out = run_filt(xyz_injected.traj[ind], xyz_injected.ins[ind], 
                        getfield(xyz_injected, cfg.mag_use)[ind], itp_mapS,
                        :ekf;
                        R    = cfg.R)

    # -- Extract position errors -----------------------------------------------
    # run_filt returns FILTout; lat/lon errors are in metres
    lat_est   = filt_out.lat
    lon_est   = filt_out.lon
    lat_true  = xyz_orig.traj.lat[ind]
    lon_true  = xyz_orig.traj.lon[ind]

    err_north, err_east = latlon_to_ne_error(lat_est, lon_est, lat_true, lon_true)

    t       = xyz_orig.traj.tt[ind]
    heading = xyz_orig.ins_yaw[ind]
    comp_mag = getfield(xyz_injected, cfg.mag_use)[ind]

    verbose && @printf "  ✓  DRMS = %.1f m\n" sqrt(mean(err_north.^2 + err_east.^2))

    return BenchmarkResult(
        label,
        Float64.(t .- t[1]),   # normalise to seconds from segment start
        lat_est, lon_est,
        lat_true, lon_true,
        err_north, err_east,
        Float64.(heading),
        Float64.(comp_mag),
        filt_out,
        cfg
    )
end

# -- Coordinate conversion helper ---------------------------------------------

"""
    latlon_to_ne_error(lat_est, lon_est, lat_true, lon_true) -> (north_m, east_m)

Convert lat/lon position errors to northing/easting in metres using a
local tangent-plane approximation (accurate to < 0.1 m over 100 km at mid-lat).

Earth radius: WGS-84 mean radius = 6,371,000 m.
"""
function latlon_to_ne_error(lat_est, lon_est, lat_true, lon_true)
    R_earth = 6_371_000.0        # mean radius (m)
    lat_mid = mean(lat_true)

    Δlat = lat_est .- lat_true
    Δlon = lon_est .- lon_true

    north_m = R_earth .* Δlat
    east_m  = R_earth .* cos(lat_mid) .* Δlon

    return north_m, east_m
end

# -- Convenience: run both at once ---------------------------------------------

"""
    run_both_benchmarks(data, cal_ind, eval_ind, custom_model, mapS;
                        config=default_config(), verbose=true)
                        -> (result_official, result_custom)

High-level entry point: compensates both official and custom TL, injects them,
and returns two `BenchmarkResult`s ready for metric computation and plotting.

# Arguments
- `data`         : `BenchmarkData` struct (from `load_benchmark_data`)
- `cal_ind`      : BitVector — calibration segment
- `eval_ind`     : BitVector — evaluation segment (must not overlap cal_ind)
- `custom_model` : fitted `CustomTLModel`
- `mapS`         : MagNav MapS
"""
function run_both_benchmarks(
    data         :: Any,     # BenchmarkData
    cal_ind      :: AbstractVector{Bool},
    eval_ind     :: AbstractVector{Bool},
    custom_model :: CustomTLModel,
    mapS         :: Any;
    config       :: Any    = nothing,
    verbose      :: Bool   = true
)
    cfg = isnothing(config) ? default_config() : config
    xyz = data.xyz

    # -- Validate segment non-overlap -----------------------------------------
    overlap = sum(cal_ind .& eval_ind)
    overlap > 0 && @warn "cal_ind and eval_ind overlap by $overlap samples — " *
                         "this will inflate performance metrics."

    # -- Official TL -----------------------------------------------------------
    verbose && println("\n-- Preparing Official TL compensation --")
    verbose && println("  DEBUG: cfg.tl_terms = ", cfg.tl_terms, " (type: ", typeof(cfg.tl_terms), ")")
    official = compensate_full_flight_official(
        xyz, cal_ind, eval_ind;
        terms     = cfg.tl_terms,
        mag_field = cfg.mag_use,
        flux_field= cfg.flux_use
    )
    xyz_official = inject_compensation(xyz, official.signal; mag_field=cfg.mag_use)
    result_official = run_ekf_benchmark(
        "official_tl", xyz_official, xyz, mapS, eval_ind;
        config=cfg, verbose=verbose
    )

    # -- Custom TL -------------------------------------------------------------
    verbose && println("\n-- Preparing Custom TL compensation --")
    comp_custom_full = compensate_full_flight(
        custom_model, xyz, cal_ind, eval_ind;
        detrend_window = cfg.detrend_window,
        mag_field      = cfg.mag_use,
        flux_field     = cfg.flux_use
    )
    xyz_custom = inject_compensation(xyz, comp_custom_full; mag_field=cfg.mag_use)
    result_custom = run_ekf_benchmark(
        "custom_tl", xyz_custom, xyz, mapS, eval_ind;
        config=cfg, verbose=verbose
    )

    return (official=result_official, custom=result_custom)
end
