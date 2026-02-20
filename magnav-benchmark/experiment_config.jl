"""
    experiment_config.jl

Typed configuration for all tunable benchmark parameters.

Design principle: every numeric constant that affects results must live here.
The config is serialised alongside every result set so that any run can be
exactly reproduced.
"""

# -- Config type ---------------------------------------------------------------

"""
    ExperimentConfig

All parameters that can affect benchmark results.

EKF noise parameters follow MagNav.jl conventions:
- P0 : initial state covariance (diagonal, position in rad²)
- Qd : process noise covariance
- R  : measurement noise variance (nT²)
"""
mutable struct ExperimentConfig
    # -- Sensor assignments -------------------------------------------------
    mag_use    :: Symbol     # scalar mag field, e.g. :mag_1_c
    flux_use   :: Symbol     # fluxgate field, e.g. :flux_a

    # -- TL model ----------------------------------------------------------
    tl_terms   :: Symbol     # MagNav term set for official TL
    tl_lambda  :: Float64    # ridge regularization for custom TL (0 = OLS)
    detrend_window :: Symbol # :cal, :local, or :none

    # -- EKF noise (MagNav defaults are reasonable starting points) ---------
    P0         :: Matrix{Float64}
    Qd         :: Matrix{Float64}
    R          :: Float64    # nT²

    # -- Experiment control -------------------------------------------------
    rng_seed   :: Int        # for any stochastic element
    results_dir:: String     # output directory for this run
    run_id     :: String     # unique identifier for this experiment run
    notes      :: String     # free text
end

# -- Constructor ---------------------------------------------------------------

"""
    default_config(; results_dir="results", notes="") -> ExperimentConfig

Return a sensible default configuration.

EKF state vector assumed: [δlat, δlon, δalt, δvN, δvE, δalt_bias, ...]
Adjust P0/Qd/R to match your MagNav XYZ struct variant (XYZ0/XYZ1/XYZ20).
"""
function default_config(;
    results_dir:: String = "results",
    notes      :: String = ""
) :: ExperimentConfig

    # These match MagNav.jl's own demo defaults for SGL 2020
    # Tune conservatively — larger R makes EKF more reliant on INS
    P0_diag = [1e-6, 1e-6, 1e2,    # position (rad², rad², m²)
               1e-2, 1e-2, 1e-2,   # velocity (m²/s²)
               1e4]                 # mag bias (nT²)
    Qd_diag = [1e-10, 1e-10, 1e-1,
               1e-4,  1e-4,  1e-4,
               1e1]
    P0 = diagm(P0_diag)
    Qd = diagm(Qd_diag)

    run_id = Dates.format(now(), "yyyymmdd_HHMMSS")

    return ExperimentConfig(
        :mag_1_c,                  # mag_use
        :flux_a,                   # flux_use
        :permanent_plus_induced,   # tl_terms
        0.0,                       # tl_lambda (OLS)
        :local,                    # detrend_window
        P0,                        # EKF initial cov
        Qd,                        # EKF process noise
        100.0,                     # R (nT²) — measurement noise
        42,                        # rng_seed
        results_dir,
        run_id,
        notes
    )
end

# -- Persistence ---------------------------------------------------------------

"""
    save_config(cfg::ExperimentConfig, path::String)

Save configuration to a human-readable TOML-like text file.
(No external TOML dep — uses plain write for portability.)
"""
function save_config(cfg::ExperimentConfig, path::String)
    isdir(dirname(path)) || mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# ExperimentConfig — generated $(now())")
        println(io, "run_id          = \"$(cfg.run_id)\"")
        println(io, "mag_use         = :$(cfg.mag_use)")
        println(io, "flux_use        = :$(cfg.flux_use)")
        println(io, "tl_terms        = :$(cfg.tl_terms)")
        println(io, "tl_lambda       = $(cfg.tl_lambda)")
        println(io, "detrend_window  = :$(cfg.detrend_window)")
        println(io, "R               = $(cfg.R)")
        println(io, "rng_seed        = $(cfg.rng_seed)")
        println(io, "results_dir     = \"$(cfg.results_dir)\"")
        println(io, "notes           = \"$(cfg.notes)\"")
        println(io, "P0_diag         = $(diag(cfg.P0))")
        println(io, "Qd_diag         = $(diag(cfg.Qd))")
    end
    println("  Config saved: $path")
end

"""
    load_config(path::String) -> ExperimentConfig

Reload a config from a saved text file.
(Parses key = value lines; complex matrices fall back to defaults.)
"""
function load_config(path::String)
    isfile(path) || error("Config file not found: $path")
    cfg = default_config()

    for line in eachline(path)
        s = strip(line)
        (startswith(s, "#") || isempty(s)) && continue
        parts = split(line, "="; limit=2)
        length(parts) < 2 && continue
        key = strip(parts[1])
        val = strip(parts[2])

        key == "run_id"         && (cfg.run_id         = String(strip(val, '"')))
        key == "tl_lambda"      && (cfg.tl_lambda      = parse(Float64, val))
        key == "R"              && (cfg.R               = parse(Float64, val))
        key == "rng_seed"       && (cfg.rng_seed        = parse(Int, val))
        key == "notes"          && (cfg.notes           = String(strip(val, '"')))
        key == "results_dir"    && (cfg.results_dir     = String(strip(val, '"')))
    end

    return cfg
end

# -- Utility -------------------------------------------------------------------

"""
    describe_config(cfg::ExperimentConfig)

Pretty-print the current configuration.
"""
function describe_config(cfg::ExperimentConfig)
    @printf "\n-- ExperimentConfig ---------------------------------\n"
    @printf "  Run ID         : %s\n"  cfg.run_id
    @printf "  mag_use        : %s\n"  string(cfg.mag_use)
    @printf "  flux_use       : %s\n"  string(cfg.flux_use)
    @printf "  tl_terms       : %s\n"  string(cfg.tl_terms)
    @printf "  tl_lambda      : %.2e\n" cfg.tl_lambda
    @printf "  detrend_window : %s\n"  string(cfg.detrend_window)
    @printf "  R (EKF)        : %.1f nT²\n" cfg.R
    @printf "  rng_seed       : %d\n"  cfg.rng_seed
    @printf "  results_dir    : %s\n"  cfg.results_dir
    @printf "  notes          : %s\n"  (isempty(cfg.notes) ? "(none)" : cfg.notes)
    @printf "----------------------------------------------------\n\n"
end
