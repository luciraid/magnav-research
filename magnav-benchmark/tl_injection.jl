"""
    tl_injection.jl

Injection layer: bridges a compensated scalar signal into MagNav's EKF pipeline
without modifying MagNav source code.

Strategy
--------
MagNav's `run_filt` (specifically the EKF variant) accepts an XYZ struct whose
`mag_1_c` field contains the scalar magnetometer readings.  We create a *shallow
copy* of the XYZ struct with only `mag_1_c` (and optionally `mag_1_uc`) replaced
by the custom-compensated signal.  Everything else — INS, GPS truth, map, etc. —
remains identical to the original object.

This guarantees:
  1. The INS trajectory is identical for both runs (no navigation drift difference)
  2. The map interpolant is identical
  3. The only independent variable is the compensation signal

Copy strategy
-------------
MagNav XYZ structs (XYZ0, XYZ1, XYZ20) are plain Julia structs with named fields.
We use `Base.setproperty` idiom via `@set` (Accessors.jl) if available, falling
back to a manual copy constructor that works with field reflection.
"""

# -- Shallow XYZ copy with replaced mag field ----------------------------------

"""
    inject_compensation(xyz, compensated_full; mag_field=:mag_1_c) -> xyz_new

Return a copy of `xyz` where `mag_field` is replaced with `compensated_full`.

`compensated_full` must be the same length as the full flight (not just the
calibration segment) — it is the compensated signal for ALL indices that will
be passed to `run_filt`.

# Arguments
- `xyz`              : original MagNav XYZ struct
- `compensated_full` : full-length compensated scalar vector (nT)
- `mag_field`        : symbol of the field to replace (default `:mag_1_c`)

# Returns
- New struct of the same concrete type as `xyz`, with `mag_field` replaced.
"""
function inject_compensation(
    xyz              :: T,
    compensated_full :: AbstractVector{<:Real};
    mag_field        :: Symbol = :mag_1_c
) where T
    n_full = length(getfield(xyz, mag_field))
    n_comp = length(compensated_full)

    if n_full != n_comp
        error("""
        inject_compensation: length mismatch.
          xyz.$mag_field has $n_full samples,
          compensated_full has $n_comp samples.
        compensated_full must cover the ENTIRE flight, not just the calibration
        window.  Use `compensate_full_flight` to produce a full-length signal.
        """)
    end

    return _copy_with_field(xyz, mag_field, Float64.(compensated_full))
end

"""
    _copy_with_field(xyz, field, new_value) -> xyz_new

Construct a copy of `xyz` with one field replaced, using Julia's field
reflection.  Works for any concrete struct type that supports a constructor
from all its fields in order.

This is the core compatibility shim — it does not use Accessors.jl to avoid
adding an unresolvable dependency.
"""
function _copy_with_field(xyz::T, field::Symbol, new_value) where T
    fields    = fieldnames(T)
    new_vals  = map(fields) do f
        f == field ? new_value : getfield(xyz, f)
    end
    return T(new_vals...)
end

# -- Full-flight compensation ---------------------------------------------------

"""
    compensate_full_flight(model::CustomTLModel, xyz, cal_ind, eval_ind;
                           detrend_window=:cal) -> Vector{Float64}

Produce a full-length compensated scalar vector for the entire flight by
applying the custom TL model trained on `cal_ind` to the whole time series.

Two detrending modes:
- `:cal`   — fit the linear trend on the calibration segment, extrapolate to full
- `:local` — fit a separate trend per evaluation segment (more robust to drift)
- `:none`  — skip detrending (not recommended unless signal is pre-detrended)

# Returns
- Vector of length `length(xyz.mag_1_c)` with NaN outside both cal and eval.
  Compensated values are written only for `eval_ind`; all others are left as
  the original detrended scalar so that MagNav can skip them.
"""
function compensate_full_flight(
    model         :: CustomTLModel,
    xyz           :: Any,
    cal_ind       :: AbstractVector{Bool},
    eval_ind      :: AbstractVector{Bool};
    detrend_window:: Symbol = :cal,
    mag_field     :: Symbol = :mag_1_c,
    flux_field    :: Symbol = :flux_a
)
    @assert model.fitted "Fit the model before calling compensate_full_flight"

    mag_full  = Float64.(getfield(xyz, mag_field))
    flux      = getfield(xyz, flux_field)
    Bx_full   = Float64.(flux.x)
    By_full   = Float64.(flux.y)
    Bz_full   = Float64.(flux.z)
    n         = length(mag_full)

    # Start with original (will overwrite eval segment)
    output = copy(mag_full)

    # -- Detrend ---------------------------------------------------------------
    t = collect(1.0:n)

    if detrend_window == :cal
        # Fit trend on calibration window, apply globally
        t_cal   = t[cal_ind]
        y_cal   = mag_full[cal_ind]
        tc      = [ones(length(t_cal)) t_cal] \ y_cal
        trend   = [ones(n) t] * tc
        y_dt    = mag_full .- trend

    elseif detrend_window == :local
        # Fit trend on eval window (forward-compatible with any segment)
        t_ev  = t[eval_ind]
        y_ev  = mag_full[eval_ind]
        tc    = [ones(length(t_ev)) t_ev] \ y_ev
        trend_ev = [ones(length(t_ev)) t_ev] * tc

        y_dt          = copy(mag_full)
        y_dt[eval_ind] = mag_full[eval_ind] .- trend_ev

    elseif detrend_window == :none
        y_dt = copy(mag_full)
    else
        error("detrend_window must be :cal, :local, or :none; got :$detrend_window")
    end

    # -- Apply model to eval segment -------------------------------------------
    A_eval        = build_A9(Bx_full[eval_ind], By_full[eval_ind], Bz_full[eval_ind])
    interference  = A_eval * model.coef
    output[eval_ind] = y_dt[eval_ind] .- interference

    return output
end

# -- Analogous helper for official TL -----------------------------------------

"""
    compensate_full_flight_official(xyz, cal_ind, eval_ind;
                                    terms=:permanent_plus_induced)
                                    -> (coef, Vector{Float64})

Train official MagNav TL on `cal_ind`, apply to `eval_ind`, return a
full-length compensation vector (same format as `compensate_full_flight`).
"""
function compensate_full_flight_official(
    xyz      :: Any,
    cal_ind  :: AbstractVector{Bool},
    eval_ind :: AbstractVector{Bool};
    terms    :: Symbol = :permanent_plus_induced,
    mag_field:: Symbol = :mag_1_c,
    flux_field:: Symbol = :flux_a
)
    println("DEBUG compensate_full_flight_official: terms=$(terms), type=$(typeof(terms))")
    flux    = getfield(xyz, flux_field)
    mag_full = Float64.(getfield(xyz, mag_field))
    n        = length(mag_full)

    # Fit on calibration window
    # Map high-level tl terms (symbols) to the vector-of-symbols expected
    # by MagNav.create_TL_A
    terms_arg = isa(terms, Symbol) ? (
        terms == :permanent_plus_induced ? [:permanent, :induced] : [terms]
    ) : terms
    println("DEBUG: About to call create_TL_A with terms=$(terms_arg), type=$(typeof(terms_arg))")
    A_cal   = create_TL_A(flux.x[cal_ind], flux.y[cal_ind], flux.z[cal_ind];
                           terms=terms_arg)
    mag_dt_cal = detrend(mag_full[cal_ind])
    # MagNav's create_TL_coef signature expects (Bx,By,Bz,B) or MagV; the
    # benchmark previously called a version that accepted (A, y). Implement a
    # lightweight fallback: solve least-squares (or ridge when λ > 0).
    λ = 0.0
    if λ == 0.0
        coef = A_cal \ mag_dt_cal
    else
        coef = (A_cal' * A_cal + λ * I) \ (A_cal' * mag_dt_cal)
    end

    # Apply to eval window (detrend locally)
    t         = collect(1.0:n)
    t_ev      = t[eval_ind]
    y_ev      = mag_full[eval_ind]
    tc        = [ones(length(t_ev)) t_ev] \ y_ev
    trend_ev  = [ones(length(t_ev)) t_ev] * tc
    y_dt_ev   = y_ev .- trend_ev

    A_eval   = create_TL_A(flux.x[eval_ind], flux.y[eval_ind], flux.z[eval_ind];
                            terms=terms_arg)
    comp_ev  = y_dt_ev .- A_eval * coef

    output           = copy(mag_full)
    output[eval_ind] = comp_ev

    return (coef=coef, signal=output)
end
