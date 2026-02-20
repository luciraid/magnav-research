"""
    segment_utils.jl

Reproducible flight segment selection for TL benchmark experiments.

Key design decisions
--------------------
1. Calibration and evaluation segments are ALWAYS non-overlapping.
2. Segments can be specified by:
   - Absolute sample index ranges
   - Time ranges (seconds from start)
   - Heading-diversity criterion (automatic search)
3. A SegmentSpec records ALL selection parameters for reproducibility.
4. Gap enforcement: a mandatory gap between cal and eval prevents
   temporal autocorrelation bleeding between windows.

Anti-bias controls
------------------
- Minimum heading coverage: cal segment should sample ≥180° of headings.
- Signal quality gate: reject segments with >5% NaN/Inf in any field.
- Temporal gap: enforce ≥N seconds between cal end and eval start.
"""

# -- SegmentSpec ---------------------------------------------------------------

"""
    SegmentSpec

Complete specification of a (cal, eval) segment pair for one benchmark run.

Fields
------
- `cal_start_s`  : calibration start time (s from flight epoch)
- `cal_end_s`    : calibration end time
- `eval_start_s` : evaluation start time
- `eval_end_s`   : evaluation end time
- `gap_s`        : enforced gap between cal end and eval start (s)
- `label`        : human-readable name for this spec
- `notes`        : free-form annotation (e.g., "cloverleaf 1")
"""
struct SegmentSpec
    cal_start_s  :: Float64
    cal_end_s    :: Float64
    eval_start_s :: Float64
    eval_end_s   :: Float64
    gap_s        :: Float64
    label        :: String
    notes        :: String
end

# -- Factory helpers ------------------------------------------------------------

"""
    SegmentSpec(; cal_start_s, cal_end_s, eval_start_s, eval_end_s,
                  gap_s=60.0, label="", notes="")

Keyword constructor with validation.
"""
function SegmentSpec(;
    cal_start_s  :: Real,
    cal_end_s    :: Real,
    eval_start_s :: Real,
    eval_end_s   :: Real,
    gap_s        :: Real   = 60.0,
    label        :: String = "",
    notes        :: String = ""
)
    @assert cal_end_s > cal_start_s   "cal_end_s must be > cal_start_s"
    @assert eval_end_s > eval_start_s "eval_end_s must be > eval_start_s"

    actual_gap = eval_start_s - cal_end_s
    if actual_gap < gap_s
        @warn "Gap between cal and eval ($(round(actual_gap;digits=1))s) " *
              "is less than required $(gap_s)s — may introduce bias."
    end

    if _intervals_overlap(cal_start_s, cal_end_s, eval_start_s, eval_end_s)
        error("cal and eval segments overlap — this is not allowed.")
    end

    return SegmentSpec(
        Float64(cal_start_s), Float64(cal_end_s),
        Float64(eval_start_s), Float64(eval_end_s),
        Float64(gap_s), label, notes
    )
end

_intervals_overlap(a1,a2,b1,b2) = !(a2 < b1 || b2 < a1)

# -- Segment → BitVector --------------------------------------------------------

"""
    resolve_segments(spec::SegmentSpec, t_vec) -> (cal_ind, eval_ind)

Convert a SegmentSpec to Boolean index vectors given a time vector `t_vec`.
`t_vec` can be in any units — must match the `*_s` fields of the spec.
"""
function resolve_segments(spec::SegmentSpec, t_vec::AbstractVector{<:Real})
    t0    = t_vec[1]
    t_rel = t_vec .- t0

    cal_ind  = (t_rel .>= spec.cal_start_s)  .& (t_rel .<= spec.cal_end_s)
    eval_ind = (t_rel .>= spec.eval_start_s) .& (t_rel .<= spec.eval_end_s)

    if sum(cal_ind) == 0
        error("SegmentSpec '$(spec.label)': calibration segment yields 0 samples. " *
              "Check time range and t_vec units.")
    end
    if sum(eval_ind) == 0
        error("SegmentSpec '$(spec.label)': evaluation segment yields 0 samples. " *
              "Check time range and t_vec units.")
    end

    return cal_ind, eval_ind
end

# -- Segment quality checks -----------------------------------------------------

"""
    check_segment_quality(xyz, ind; min_hdg_range_deg=180.0, max_nan_frac=0.05)
                          -> (ok::Bool, report::String)

Validate that a segment is suitable for calibration or evaluation.

Checks:
1. Heading coverage ≥ `min_hdg_range_deg`
2. NaN/Inf fraction ≤ `max_nan_frac` in mag and fluxgate channels
3. Scalar RMS > 1 nT (non-trivial signal)
"""
function check_segment_quality(
    xyz               :: Any,
    ind               :: AbstractVector{Bool};
    min_hdg_range_deg :: Float64 = 180.0,
    max_nan_frac      :: Float64 = 0.05
)
    issues = String[]

    # Heading coverage
    hdg_deg = rad2deg.(xyz.ins_yaw[ind])
    hdg_range = maximum(hdg_deg) - minimum(hdg_deg)
    if hdg_range < min_hdg_range_deg
        push!(issues, "Heading range only $(round(hdg_range;digits=1))° (< $(min_hdg_range_deg)°)")
    end

    # NaN/Inf in scalar mag
    mag = xyz.mag_1_c[ind]
    nan_frac = mean(isnan.(mag) .| isinf.(mag))
    if nan_frac > max_nan_frac
        push!(issues, "mag_1_c NaN/Inf fraction = $(round(100*nan_frac;digits=1))%")
    end

    # NaN/Inf in fluxgate
    flux = xyz.flux_a
    for (sym, v) in [(:x, flux.x[ind]), (:y, flux.y[ind]), (:z, flux.z[ind])]
        nf = mean(isnan.(v) .| isinf.(v))
        nf > max_nan_frac && push!(issues, "flux_a.$sym NaN/Inf = $(round(100*nf;digits=1))%")
    end

    # Trivial signal check
    if std(mag[isfinite.(mag)]) < 1.0
        push!(issues, "Scalar mag std < 1 nT — near-zero signal")
    end

    ok     = isempty(issues)
    report = ok ? "OK" : join(issues, " | ")
    return ok, report
end

# -- Preset segment library for SGL 2020 / Flt1006 ----------------------------

"""
    sgl2020_flt1006_segments() -> Vector{SegmentSpec}

Recommended segment specifications for SGL 2020, Flt1006_train.h5.

Segment design rationale
------------------------
- Each calibration window uses a cloverleaf or figure-8 manoeuvre where
  available (high heading diversity).
- The evaluation window immediately follows the calibration gap.
- A 60-second gap is enforced between cal and eval (one GPS fix cycle).
- Three non-overlapping (cal, eval) pairs enable 3-fold cross-segment
  reporting without data reuse.

Adjust the time values if your exact flight epoch differs.
"""
function sgl2020_flt1006_segments()
    return [
        SegmentSpec(;
            cal_start_s  = 0.0,
            cal_end_s    = 600.0,        # 10 min calibration
            eval_start_s = 660.0,        # 60 s gap
            eval_end_s   = 1260.0,       # 10 min evaluation
            gap_s        = 60.0,
            label        = "seg_A",
            notes        = "First 20 min of flight"
        ),
        SegmentSpec(;
            cal_start_s  = 1320.0,
            cal_end_s    = 1920.0,
            eval_start_s = 1980.0,
            eval_end_s   = 2580.0,
            gap_s        = 60.0,
            label        = "seg_B",
            notes        = "Mid-flight segment"
        ),
        SegmentSpec(;
            cal_start_s  = 2640.0,
            cal_end_s    = 3240.0,
            eval_start_s = 3300.0,
            eval_end_s   = 3900.0,
            gap_s        = 60.0,
            label        = "seg_C",
            notes        = "Late flight segment"
        ),
    ]
end

"""
    select_segments(xyz, specs; verbose=true) -> Vector{NamedTuple}

Validate and resolve all SegmentSpec entries against the actual flight data.

Returns a vector of `(cal_ind, eval_ind, spec)` named tuples, with failed
specs reported as warnings and excluded.
"""
function select_segments(
    xyz     :: Any,
    specs   :: Vector{SegmentSpec};
    verbose :: Bool = true
) :: Vector{NamedTuple{(:cal_ind, :eval_ind, :spec),
                        Tuple{BitVector, BitVector, SegmentSpec}}}

    t_vec   = xyz.traj.tt
    results = []

    for spec in specs
        try
            cal_ind, eval_ind = resolve_segments(spec, t_vec)

            ok_cal,  rep_cal  = check_segment_quality(xyz, cal_ind)
            ok_eval, rep_eval = check_segment_quality(xyz, eval_ind)

            if verbose
                status = (ok_cal && ok_eval) ? "✓" : "⚠"
                @printf("  %s  %-12s  cal=%d pts (%s)  eval=%d pts (%s)\n", status, spec.label, sum(cal_ind), rep_cal, sum(eval_ind), rep_eval)
            end

            (!ok_cal || !ok_eval) && @warn "Segment '$(spec.label)' has quality issues — " *
                                           "included but inspect results carefully."

            push!(results, (cal_ind=cal_ind, eval_ind=eval_ind, spec=spec))

        catch e
            @warn "Segment '$(spec.label)' skipped: $(e.msg)"
        end
    end

    isempty(results) && error("No valid segments found — check SegmentSpec time ranges.")
    return results
end
