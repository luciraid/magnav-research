### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 8b5b5b5b-5b5b-4b5b-8b5b-5b5b5b5b5b5b
begin
	# MagNav Benchmark Framework - Pluto Edition
	# Comparing Custom 9-term vs Official Tolles-Lawson Compensation

	using Pkg
	Pkg.activate(".")

	# Core dependencies
	using MagNav
	using DataFrames
	using LinearAlgebra
	using Statistics
	using Printf
	using Plots
	using PlutoUI
end

# ╔═╡ 5b5b5b5b-5b5b-4b5b-8b5b-5b5b5b5b5b5b
md"""
# MagNav TL Benchmark Framework

**Interactive Pluto Notebook**

This notebook compares custom 9-term Tolles-Lawson (TL) magnetic compensation against MagNav's official implementation using SGL 2020 flight data and Extended Kalman Filter (EKF) navigation performance evaluation.

## Configuration
"""

# ╔═╡ 4b5b5b5b-5b5b-4b5b-8b5b-5b5b5b5b5b5b
begin
	# Experiment Configuration
	@kwdef mutable struct ExperimentConfig
	    mag_use    :: Symbol     = :mag_1_c
	    flux_use   :: Symbol     = :flux_a
	    tl_terms   :: Symbol     = :permanent_plus_induced
	    tl_lambda  :: Float64    = 0.0
	    detrend_window :: Symbol = :local
	    R          :: Float64    = 100.0
	    rng_seed   :: Int        = 42
	end

	cfg = ExperimentConfig()
end

# ╔═╡ 3b5b5b5b-5b5b-4b5b-8b5b-5b5b5b5b5b5b
md"""
## Data Loading

Load SGL 2020 flight data and magnetic anomaly maps.
"""

# ╔═╡ 2b5b5b5b-5b5b-4b5b-8b5b-5b5b5b5b5b5b
begin
	# Load flight data
	flight = "Flt1006"
	xyz = get_XYZ20(flight; silent=true)
	mapS = get_map(MagNav.ottawa_area_maps; silent=true)

	md"""
	**Flight Data Loaded:**
	- Flight: $(flight)
	- Samples: $(length(xyz.traj.lat))
	- Duration: $(round(length(xyz.traj.lat)/10/60, digits=1)) minutes
	- Lat range: $(round(minimum(xyz.traj.lat)*180/π, digits=4))° to $(round(maximum(xyz.traj.lat)*180/π, digits=4))°
	- Lon range: $(round(minimum(xyz.traj.lon)*180/π, digits=4))° to $(round(maximum(xyz.traj.lon)*180/π, digits=4))°
	"""
end

# ╔═╡ 1b5b5b5b-5b5b-4b5b-8b5b-5b5b5b5b5b5b
md"""
## Segment Selection

Select calibration and evaluation segments for TL fitting and EKF testing.
"""

# ╔═╡ 0b5b5b5b-5b5b-4b5b-8b5b-5b5b5b5b5b5b
begin
	# Segment definitions (simplified for notebook)
	segments = Dict(
	    "seg_A" => (cal_start=1, cal_end=6001, eval_start=6002, eval_end=10622),
	    "seg_C" => (cal_start=20000, cal_end=26001, eval_start=26002, eval_end=32002)
	)

	selected_segment = "seg_A"
	seg = segments[selected_segment]

	cal_ind = seg.cal_start:seg.cal_end
	eval_ind = seg.eval_start:seg.eval_end

	md"""
	**Selected Segment: $(selected_segment)**
	- Calibration: $(length(cal_ind)) samples
	- Evaluation: $(length(eval_ind)) samples
	"""
end

# ╔═╡ 9b5b5b5b-5b5b-4b5b-8b5b-5b5b5b5b5b5b
md"""
## Custom TL Model Fitting

Fit the 9-term Tolles-Lawson compensation model.
"""

# ╔═╡ ab5b5b5b-5b5b-4b5b-8b5b-5b5b5b5b5b5b
begin
	# Custom 9-term TL model
	function fit_custom_TL(xyz::MagNav.XYZ20, ind::UnitRange)
	    # Extract calibration data
	    lat = xyz.traj.lat[ind]
	    lon = xyz.traj.lon[ind]
	    alt = xyz.traj.alt[ind]
	    yaw = xyz.ins_yaw[ind]
	    pitch = xyz.ins_pitch[ind]
	    roll = xyz.ins_roll[ind]
	    mag = getfield(xyz, cfg.mag_use)[ind]
	    flux = getfield(xyz, cfg.flux_use)[ind]

	    # Design matrix for 9-term model
	    A = hcat(
	        ones(length(ind)),      # constant
	        cos.(yaw), sin.(yaw),   # permanent
	        cos.(pitch), sin.(pitch), # induced
	        cos.(roll), sin.(roll)    # eddy
	    )

	    # Fit model
	    coef = A \ mag
	    mag_pred = A * coef
	    residual = mag - mag_pred

	    # Statistics
	    rms_resid = sqrt(mean(residual.^2))
	    var_reduction = (1 - var(residual)/var(mag)) * 100
	    corr_resid_hdg = cor(residual, yaw)

	    return (coef=coef, rms=rms_resid, var_reduction=var_reduction,
	            corr=corr_resid_hdg, residual=residual, pred=mag_pred)
	end

	custom_model = fit_custom_TL(xyz, cal_ind)
end

# ╔═╡ bb5b5b5b-5b5b-4b5b-8b5b-5b5b5b5b5b5b
md"""
**Custom TL Fit Results:**
- RMS Residual: $(round(custom_model.rms, digits=2)) nT
- Variance Reduction: $(round(custom_model.var_reduction, digits=1))%
- Correlation with Heading: $(round(custom_model.corr, digits=3))
"""

# ╔═╡ cb5b5b5b-5b5b-4b5b-8b5b-5b5b5b5b5b5b
begin
	# Plot TL fit results
	p1 = plot(xyz.traj.lat[cal_ind]*180/π, getfield(xyz, cfg.mag_use)[cal_ind],
	          label="Measured", linewidth=1, alpha=0.7)
	plot!(p1, xyz.traj.lat[cal_ind]*180/π, custom_model.pred,
	      label="Custom TL Fit", linewidth=2, color=:red)

	p2 = plot(xyz.traj.lat[cal_ind]*180/π, custom_model.residual,
	          label="Residual", linewidth=1, color=:blue)
	hline!(p2, [0], linestyle=:dash, color=:black, label="")

	plot(p1, p2, layout=(2,1), size=(800,400),
	     xlabel="Latitude [°]", ylabel=["Magnetic Field [nT]", "Residual [nT]"])
end

# ╔═╡ db5b5b5b-5b5b-4b5b-8b5b-5b5b5b5b5b5b
md"""
## Official TL Compensation

Apply MagNav's official TL compensation for comparison.
"""

# ╔═╡ eb5b5b5b-5b5b-4b5b-8b5b-5b5b5b5b5b5b
begin
	# Official TL compensation
	function apply_official_TL(xyz::MagNav.XYZ20, terms::Symbol)
	    # Convert symbol to vector for create_TL_A
	    term_dict = Dict(
	        :permanent => [:permanent],
	        :induced => [:induced],
	        :permanent_plus_induced => [:permanent, :induced],
	        :full => [:permanent, :induced, :eddy, :bias]
	    )
	    terms_vec = get(term_dict, terms, [:permanent, :induced])

	    # Create design matrix
	    A = create_TL_A(xyz.traj.lat, xyz.traj.lon, xyz.traj.alt,
	                    xyz.ins_yaw, xyz.ins_pitch, xyz.ins_roll;
	                    terms=terms_vec)

	    # Fit coefficients
	    mag = getfield(xyz, cfg.mag_use)
	    coef = A \ mag

	    # Apply compensation
	    mag_comp = mag - A * coef

	    return mag_comp, coef
	end

	mag_official, coef_official = apply_official_TL(xyz, cfg.tl_terms)
end

# ╔═╡ fb5b5b5b-5b5b-4b5b-8b5b-5b5b5b5b5b5b
md"""
## TL Compensation Comparison

Compare the performance of custom vs official TL compensation.
"""

# ╔═╡ gb5b5b5b-5b5b-4b5b-8b5b-5b5b5b5b5b5b
begin
	# Apply custom compensation to full flight
	function apply_custom_compensation(xyz::MagNav.XYZ20, model)
	    yaw = xyz.ins_yaw
	    pitch = xyz.ins_pitch
	    roll = xyz.ins_roll
	    mag = getfield(xyz, cfg.mag_use)

	    A_full = hcat(
	        ones(length(mag)),
	        cos.(yaw), sin.(yaw),
	        cos.(pitch), sin.(pitch),
	        cos.(roll), sin.(roll)
	    )

	    mag_comp = mag - A_full * model.coef
	    return mag_comp
	end

	mag_custom = apply_custom_compensation(xyz, custom_model)

	# Compare on evaluation segment
	eval_mag_orig = getfield(xyz, cfg.mag_use)[eval_ind]
	eval_mag_custom = mag_custom[eval_ind]
	eval_mag_official = mag_official[eval_ind]

	rms_custom = sqrt(mean((eval_mag_custom .- mean(eval_mag_custom)).^2))
	rms_official = sqrt(mean((eval_mag_official .- mean(eval_mag_official)).^2))

	md"""
	**Evaluation Segment RMS (detrended):**
	- Original: $(round(sqrt(mean((eval_mag_orig .- mean(eval_mag_orig)).^2)), digits=2)) nT
	- Custom TL: $(round(rms_custom, digits=2)) nT
	- Official TL: $(round(rms_official, digits=2)) nT
	"""
end

# ╔═╡ hb5b5b5b-5b5b-4b5b-8b5b-5b5b5b5b5b5b
begin
	# Plot compensation comparison
	p3 = plot(xyz.traj.lat[eval_ind]*180/π, eval_mag_orig,
	          label="Original", linewidth=1, alpha=0.7, color=:black)
	plot!(p3, xyz.traj.lat[eval_ind]*180/π, eval_mag_custom,
	      label="Custom TL", linewidth=2, color=:blue)
	plot!(p3, xyz.traj.lat[eval_ind]*180/π, eval_mag_official,
	      label="Official TL", linewidth=2, color=:red)

	plot(p3, size=(800,300), xlabel="Latitude [°]", ylabel="Magnetic Field [nT]",
	     title="TL Compensation Comparison - Evaluation Segment")
end

# ╔═╡ ib5b5b5b-5b5b-4b5b-8b5b-5b5b5b5b5b5b
md"""
## EKF Navigation Performance

Run Extended Kalman Filter to evaluate navigation performance with different compensation methods.
"""

# ╔═╡ jb5b5b5b-5b5b-4b5b-8b5b-5b5b5b5b5b5b
begin
	# Simplified EKF wrapper for notebook
	function run_simple_ekf(xyz::MagNav.XYZ20, mag_comp::Vector, ind::UnitRange, mapS)
	    # Create XYZ with compensated magnetometer
	    xyz_comp = deepcopy(xyz)
	    setfield!(xyz_comp, cfg.mag_use, mag_comp)

	    # Build map interpolant
	    map_vals, itp_mapS = get_map_val(mapS,
	                                     xyz_comp.traj.lat[ind],
	                                     xyz_comp.traj.lon[ind],
	                                     xyz_comp.traj.alt[ind];
	                                     return_itp=true)

	    # Run EKF (simplified - may need adjustment)
	    try
	        filt_out = run_filt(xyz_comp.traj[ind], xyz_comp.ins[ind],
	                            getfield(xyz_comp, cfg.mag_use)[ind], itp_mapS,
	                            :ekf; R=cfg.R)
	        return filt_out
	    catch e
	        @warn "EKF failed: $e"
	        return nothing
	    end
	end

	# Run EKF for both methods
	ekf_custom = run_simple_ekf(xyz, mag_custom, eval_ind, mapS)
	ekf_official = run_simple_ekf(xyz, mag_official, eval_ind, mapS)
end

# ╔═╡ kb5b5b5b-5b5b-4b5b-8b5b-5b5b5b5b5b5b
md"""
## Results Summary

**TL Compensation Performance:**
- Custom 9-term model shows $(round(custom_model.var_reduction, digits=1))% variance reduction in calibration
- Evaluation segment RMS: Custom $(round(rms_custom, digits=2)) nT vs Official $(round(rms_official, digits=2)) nT

**EKF Results:**
$(ekf_custom === nothing ? "- Custom TL EKF: Failed" : "- Custom TL EKF: Completed")
$(ekf_official === nothing ? "- Official TL EKF: Failed" : "- Official TL EKF: Completed")

This Pluto notebook provides an interactive environment for exploring and comparing different magnetic compensation approaches. Modify the configuration parameters above to test different settings!
"""

# ╔═╡ 00000000-0000-0000-0000-000000000000
md"""
## Next Steps

1. **Fix EKF Implementation:** Resolve the trajectory indexing issues for proper EKF execution
2. **Add More Visualizations:** Position error plots, trajectory comparisons, statistical analysis
3. **Parameter Studies:** Test different TL terms, regularization, and EKF settings
4. **Multi-Segment Analysis:** Extend to compare performance across different flight segments
"""

# ╔═╡ Cell order:
# ╠═8b5b5b5b-5b5b-4b5b-8b5b-5b5b5b5b5b5b
# ╟─5b5b5b5b-5b5b-4b5b-8b5b-5b5b5b5b5b5b
# ╠═4b5b5b5b-5b5b-4b5b-8b5b-5b5b5b5b5b5b
# ╟─3b5b5b5b-5b5b-4b5b-8b5b-5b5b5b5b5b5b
# ╠═2b5b5b5b-5b5b-4b5b-8b5b-5b5b5b5b5b5b
# ╟─1b5b5b5b-5b5b-4b5b-8b5b-5b5b5b5b5b5b
# ╠═0b5b5b5b-5b5b-4b5b-8b5b-5b5b5b5b5b5b
# ╟─9b5b5b5b-5b5b-4b5b-8b5b-5b5b5b5b5b5b
# ╠═ab5b5b5b-5b5b-4b5b-8b5b-5b5b5b5b5b5b
# ╟─bb5b5b5b-5b5b-4b5b-8b5b-5b5b5b5b5b5b
# ╠═cb5b5b5b-5b5b-4b5b-8b5b-5b5b5b5b5b5b
# ╟─db5b5b5b-5b5b-4b5b-8b5b-5b5b5b5b5b5b
# ╠═eb5b5b5b-5b5b-4b5b-8b5b-5b5b5b5b5b5b
# ╟─fb5b5b5b-5b5b-4b5b-8b5b-5b5b5b5b5b5b
# ╠═gb5b5b5b-5b5b-4b5b-8b5b-5b5b5b5b5b5b
# ╠═hb5b5b5b-5b5b-4b5b-8b5b-5b5b5b5b5b5b
# ╟─ib5b5b5b-5b5b-4b5b-8b5b-5b5b5b5b5b5b
# ╠═jb5b5b5b-5b5b-4b5b-8b5b-5b5b5b5b5b5b
# ╟─kb5b5b5b-5b5b-4b5b-8b5b-5b5b5b5b5b5b
# ╟─00000000-0000-0000-0000-000000000000