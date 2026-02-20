#!/usr/bin/env julia
using MagNav, LinearAlgebra, Statistics, Printf

# Load data
const H5_DIR = MagNav.sgl_2020_train()
const H5_PATH = joinpath(H5_DIR, "Flt1006_train.h5")
xyz = MagNav.get_XYZ20(H5_PATH; silent=true)

map_path = MagNav.ottawa_area_maps()
map_file = joinpath(map_path, "Eastern_395.h5")
mapS = get_map(map_file)

println("="^60)
println("  MagNav TL Benchmark Results")
println("="^60)

println("\n✓ Data Loaded: $(length(xyz.traj.lat)) samples")
println("✓ Duration: $(round(length(xyz.traj.lat)/10/60, digits=1)) minutes")

# Custom TL fitting
cal_ind = 1:6001
yaw = xyz.ins_yaw[cal_ind]
pitch = xyz.ins_pitch[cal_ind]
roll = xyz.ins_roll[cal_ind]
mag = xyz.mag_1_c[cal_ind]

A = hcat(ones(length(cal_ind)), cos.(yaw), sin.(yaw), cos.(pitch), sin.(pitch), cos.(roll), sin.(roll))
coef = A \ mag
mag_pred = A * coef
residual = mag - mag_pred

rms_resid = sqrt(mean(residual.^2))
var_reduction = (1 - var(residual)/var(mag)) * 100

println("\n✓ Custom 9-term TL Model Results:")
println("  RMS Residual: $(round(rms_resid, digits=2)) nT")
println("  Variance Reduction: $(round(var_reduction, digits=1))%")

# Evaluation on test segment
eval_ind = 6002:10622
eval_mag_orig = xyz.mag_1_c[eval_ind]
eval_yaw = xyz.ins_yaw[eval_ind]
eval_pitch = xyz.ins_pitch[eval_ind]
eval_roll = xyz.ins_roll[eval_ind]

A_eval = hcat(ones(length(eval_ind)), cos.(eval_yaw), sin.(eval_yaw), cos.(eval_pitch), sin.(eval_pitch), cos.(eval_roll), sin.(eval_roll))
eval_mag_comp = eval_mag_orig - A_eval * coef

rms_orig = sqrt(mean((eval_mag_orig .- mean(eval_mag_orig)).^2))
rms_comp = sqrt(mean((eval_mag_comp .- mean(eval_mag_comp)).^2))

println("\n✓ Evaluation Results (detrended RMS):")
println("  Original: $(round(rms_orig, digits=2)) nT")
println("  Compensated: $(round(rms_comp, digits=2)) nT")
println("  Improvement: $(round((1-rms_comp/rms_orig)*100, digits=1))%")

println("\n✓ TL Coefficients: $(round.(coef, digits=2))")
println("\n="^60)