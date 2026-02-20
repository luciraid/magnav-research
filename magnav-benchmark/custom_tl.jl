"""
    custom_tl.jl

Custom 9-term Tolles–Lawson aircraft magnetic compensation.

Model:
    B_aircraft = A_9 * c + noise

Terms (9 total, no eddy):
    Permanent (3):   cos(heading), sin(heading) components as unit vector {Bx, By, Bz}
    Induced linear (3):  Bt * {Bx/Bt, By/Bt, Bz/Bt}  (i.e., {Bx, By, Bz})
    Induced quadratic (3): {Bx^2, By^2, Bz^2} / Bt

where Bx, By, Bz are fluxgate magnetometer readings (body frame),
Bt = √(Bx²+By²+Bz²).

Reference: Tolles & Lawson (1950), Leliak (1961).
"""

using Printf

# -- Type ---------------------------------------------------------------------

"""
    CustomTLModel

Holds the fitted coefficients and diagnostics for the custom 9-term TL model.

Fields
------
- `coef::Vector{Float64}` : 9 fitted coefficients
- `fitted::Bool`          : whether model has been fit
- `residual_rms::Float64` : RMS of compensation residuals (nT)
- `var_reduction::Float64`: variance reduction fraction (0–1)
- `heading_corr_pre::Float64`  : |corr(residual, heading)| before compensation
- `heading_corr_post::Float64` : |corr(residual, heading)| after compensation
- `n_samples::Int`        : number of calibration samples used
"""
mutable struct CustomTLModel
    coef            :: Vector{Float64}
    fitted          :: Bool
    residual_rms    :: Float64
    var_reduction   :: Float64
    heading_corr_pre  :: Float64
    heading_corr_post :: Float64
    n_samples       :: Int

    CustomTLModel() = new(
        zeros(9), false,
        NaN, NaN, NaN, NaN, 0
    )
end

# -- Design matrix -------------------------------------------------------------

"""
    build_A9(Bx, By, Bz) -> Matrix{Float64}

Build the 9-term TL design matrix A (n × 9).

Columns:
 1–3  permanent:          Bx/Bt,  By/Bt,  Bz/Bt
 4–6  linear induced:     Bx,     By,     Bz
 7–9  quadratic induced:  Bx²/Bt, By²/Bt, Bz²/Bt

All columns are divided through by Bt where appropriate so that each column
has consistent units (nT or dimensionless), improving conditioning.

# Arguments
- `Bx, By, Bz` : body-frame fluxgate readings (nT), same length

# Returns
- `A9` : n × 9 design matrix
"""
function build_A9(
    Bx::AbstractVector{<:Real},
    By::AbstractVector{<:Real},
    Bz::AbstractVector{<:Real}
)
    n  = length(Bx)
    @assert length(By) == n && length(Bz) == n "Bx, By, Bz must have equal length"

    Bt = @. sqrt(Bx^2 + By^2 + Bz^2)

    # Guard against near-zero Bt (e.g., sensor dropout)
    ε = 1.0  # 1 nT floor — effectively zero in MagNav context
    Bt_safe = max.(Bt, ε)

    A9 = Matrix{Float64}(undef, n, 9)

    # Permanent (unit direction cosines)
    A9[:, 1] = Bx ./ Bt_safe
    A9[:, 2] = By ./ Bt_safe
    A9[:, 3] = Bz ./ Bt_safe

    # Linear induced (raw fluxgate)
    A9[:, 4] = Bx
    A9[:, 5] = By
    A9[:, 6] = Bz

    # Quadratic induced (scaled by 1/Bt for unit consistency)
    A9[:, 7] = Bx.^2 ./ Bt_safe
    A9[:, 8] = By.^2 ./ Bt_safe
    A9[:, 9] = Bz.^2 ./ Bt_safe

    return A9
end

# -- Fitting -------------------------------------------------------------------

"""
    fit_custom_tl!(model, Bx, By, Bz, mag_scalar, heading;
                   lambda=0.0, detrend=true, verbose=true)

Fit the custom 9-term TL model using ordinary (or ridge) least squares.

# Arguments
- `model`       : `CustomTLModel` to populate in-place
- `Bx/By/Bz`    : body-frame fluxgate components (nT)
- `mag_scalar`  : scalar magnetometer reading to be predicted (nT)
- `heading`     : aircraft true heading (radians), used for diagnostics only
- `lambda`      : ridge regularization parameter (0 = OLS)
- `detrend`     : remove linear trend from scalar before fitting (recommended)
- `verbose`     : print fit summary

# Returns
- Mutated `model` with `coef`, `fitted`, and diagnostics populated
"""
function fit_custom_tl!(
    model     :: CustomTLModel,
    Bx        :: AbstractVector{<:Real},
    By        :: AbstractVector{<:Real},
    Bz        :: AbstractVector{<:Real},
    mag_scalar:: AbstractVector{<:Real},
    heading   :: AbstractVector{<:Real};
    lambda    :: Float64  = 0.0,
    detrend   :: Bool     = true,
    verbose   :: Bool     = true
)
    n = length(mag_scalar)
    @assert length(Bx) == n "All vectors must be same length"

    # -- Pre-processing --------------------------------------------------------
    y = Float64.(mag_scalar)

    if detrend
        t = collect(1.0:n)
        trend_coef = [ones(n) t] \ y
        trend      = [ones(n) t] * trend_coef
        y          = y .- trend
    end

    # -- Build A, solve --------------------------------------------------------
    A = build_A9(Bx, By, Bz)

    if lambda ≈ 0.0
        coef = A \ y
    else
        # Ridge: (Aᵀ A + λI) coef = Aᵀ y
        coef = (A'A + lambda * I) \ (A'y)
    end

    # -- Diagnostics ----------------------------------------------------------
    residuals_post = y .- A * coef
    residuals_pre  = y             # before compensation = detrended scalar

    rms_post  = sqrt(mean(residuals_post.^2))
    var_pre   = var(residuals_pre)
    var_post  = var(residuals_post)
    var_red   = 1.0 - var_post / var_pre

    corr_pre  = abs(cor(residuals_pre,  heading))
    corr_post = abs(cor(residuals_post, heading))

    # -- Populate model --------------------------------------------------------
    model.coef             = coef
    model.fitted           = true
    model.residual_rms     = rms_post
    model.var_reduction    = var_red
    model.heading_corr_pre = corr_pre
    model.heading_corr_post= corr_post
    model.n_samples        = n

    if verbose
        @printf "\n-- Custom TL Fit Summary ---------------------------\n"
        @printf "  Samples          : %d\n"             n
        @printf "  Residual RMS     : %.4f nT\n"        rms_post
        @printf "  Variance reduction: %.1f%%\n"        100 * var_red
        @printf "  |corr(resid, hdg)|: %.4f → %.4f\n"  corr_pre corr_post
        @printf "  Coefficients     : %s\n"             string(round.(coef; digits=4))
        @printf "----------------------------------------------------\n\n"
    end

    return model
end

# -- Compensation --------------------------------------------------------------

"""
    compensate_custom_tl(model, Bx, By, Bz, mag_scalar; detrend=true) -> Vector{Float64}

Apply a fitted CustomTLModel to new data.

The output is the scalar magnetometer reading with aircraft interference removed.

# Returns
- `compensated` : vector of compensated scalar readings (nT)
"""
function compensate_custom_tl(
    model      :: CustomTLModel,
    Bx         :: AbstractVector{<:Real},
    By         :: AbstractVector{<:Real},
    Bz         :: AbstractVector{<:Real},
    mag_scalar :: AbstractVector{<:Real};
    detrend    :: Bool = true
)
    @assert model.fitted "Model must be fitted before calling compensate_custom_tl"

    n  = length(mag_scalar)
    y  = Float64.(mag_scalar)

    if detrend
        t = collect(1.0:n)
        trend_coef = [ones(n) t] \ y
        trend      = [ones(n) t] * trend_coef
        y          = y .- trend
    end

    A             = build_A9(Bx, By, Bz)
    interference  = A * model.coef
    compensated   = y .- interference

    return compensated
end

# -- Official TL helpers -------------------------------------------------------

"""
    fit_official_tl(xyz, ind; terms=:permanent_plus_induced, verbose=false)
        -> (A_matrix, coefs, compensated_scalar)

Thin wrapper around MagNav's `create_TL_A` / `create_TL_coef` / `detrend`
for use in side-by-side benchmarking.

Returns the design matrix, coefficients, and compensated scalar so the caller
can feed them into the injection pipeline.

# Arguments
- `xyz`   : MagNav XYZ data struct
- `ind`   : BitVector / index range selecting calibration segment
- `terms` : TL term set passed to `create_TL_A` (see MagNav docs)
"""
function fit_official_tl(
    xyz  :: Any,          # MagNav XYZ0 / XYZ1 / XYZ20 struct
    ind  :: AbstractVector{Bool};
    terms:: Symbol = :permanent_plus_induced,
    verbose:: Bool = false
)
    mag_use  = :mag_1_c    # strapdown scalar
    flux_use = :flux_a     # fluxgate vector sensor

    # Build design matrix
    A = create_TL_A(getfield(xyz, flux_use).x[ind],
                    getfield(xyz, flux_use).y[ind],
                    getfield(xyz, flux_use).z[ind];
                    terms = terms)

    # Detrend scalar
    mag_raw = getfield(xyz, mag_use)[ind]
    mag_dt  = detrend(mag_raw)

    # Solve for coefficients
    coef = create_TL_coef(A, mag_dt)

    if verbose
        @printf "Official TL: %d terms fitted on %d samples\n" length(coef) sum(ind)
    end

    # Compensated = detrended - predicted interference
    comp = mag_dt .- A * coef

    return (A_cal=A, coef=coef, compensated=comp)
end

"""
    compensate_official_tl(xyz, ind, coef) -> Vector{Float64}

Apply previously-fitted official TL coefficients to a (different) segment.
"""
function compensate_official_tl(
    xyz  :: Any,
    ind  :: AbstractVector{Bool},
    coef :: AbstractVector{<:Real};
    terms:: Symbol = :permanent_plus_induced
)
    flux = getfield(xyz, :flux_a)
    A    = create_TL_A(flux.x[ind], flux.y[ind], flux.z[ind]; terms=terms)
    mag  = getfield(xyz, :mag_1_c)[ind]
    mag_dt = detrend(mag)
    return mag_dt .- A * coef
end
