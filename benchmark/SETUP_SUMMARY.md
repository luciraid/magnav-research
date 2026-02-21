# MagNav Benchmark Framework - Setup Complete ✓

## Overview

This is a **Julia project for benchmarking aircraft magnetic compensation algorithms**. It compares a custom 9-term Tolles-Lawson (TL) model against the official MagNav.jl implementation at the navigation level.

## Project Structure

### Core Modules
- **custom_tl.jl** - Custom 9-term Tolles-Lawson compensation model
- **tl_injection.jl** - Inject compensation into flight data
- **ekf_wrapper.jl** - Extended Kalman Filter wrapper for benchmarking
- **nav_metrics.jl** - Navigation metrics (RMS error, DRMS, heading correlation)
- **segment_utils.jl** - Flight segment management and validation
- **data_loader.jl** - HDF5 data loading utilities
- **experiment_config.jl** - Experiment configuration

### Visualization
- **bench_plots.jl** - Matplotlib/Plots.jl visualizations

### Experiments
- **bench_tl_compare.jl** - Primary single-segment benchmark
- **bench_multi_segment.jl** - Cross-segment validation
- **bench_sensitivity.jl** - Calibration window length sweep

## How to Run

### 1. Run Tests (No Data Required)
```bash
cd /home/lucifer/magnav/magnav-research/magnav-benchmark
julia run_tests.jl
```
✓ **Status**: All 10 tests passing

### 2. Run Benchmark Experiments (Data Automatic)

**No setup required!** The HDF5 data is downloaded automatically via Julia's artifact system.

```bash
julia bench_tl_compare.jl         # Main benchmark
julia bench_multi_segment.jl      # Cross-segment validation  
julia bench_sensitivity.jl        # Parameter sensitivity
```

#### Data Location (Auto-Managed)
- On first run, `sgl_2020_train()` downloads the dataset to: `~/.julia/artifacts/<hash>/sgl_2020_train/`
- The file `Flt1006_train.h5` is automatically located there
- **You don't need to manually manage or copy files**
- Julia's artifact system keeps your repo clean and avoids committing large binaries

#### Advanced: Custom Output Directory
```bash
RESULTS_DIR=/tmp/results julia bench_tl_compare.jl
```

### Using Different Flight Files
The artifact system includes 6 flight datasets. To use a different flight (e.g., Flt1002):
1. Edit the benchmark script and change:
   ```julia
   const H5_PATH = joinpath(H5_DIR, "Flt1002_train.h5")
   ```
2. Or modify `data_loader.jl` function calls with different flight IDs

Available flights:
- Flt1002_train.h5 (149.1 MB)
- Flt1003_train.h5 (119.8 MB)
- Flt1004_train.h5 (61.1 MB)
- Flt1005_train.h5 (61.3 MB)
- **Flt1006_train.h5** (81.2 MB) ← default
- Flt1007_train.h5 (85.8 MB)

To verify artifact system is working:
```bash
julia verify_artifact_system.jl
```

## Key Classes & Functions

### CustomTLModel
```julia
struct CustomTLModel
    coef::Vector{Float64}          # 9 fitted coefficients
    fitted::Bool                   # Whether model has been fit
    residual_rms::Float64          # RMS of residuals (nT)
    var_reduction::Float64         # Variance reduction (0-1)
    heading_corr_pre::Float64      # Heading correlation before compensation
    heading_corr_post::Float64     # Heading correlation after compensation
    n_samples::Int                 # Number of calibration samples
end
```

### Main Functions
- `build_A9(Bx, By, Bz)` - Build 9-term TL design matrix
- `fit_custom_tl!(model, Bx, By, Bz, mag_scalar, heading)` - Fit the model
- `compensate_flight_data(Bx, By, Bz, coef)` - Apply compensation
- `compute_nav_metrics(Xd, X_truth)` - Compute navigation metrics

## Dependencies

All required Julia packages are installed globally:
- ✓ MagNav (official MagNav.jl)
- ✓ DataFrames
- ✓ CSV
- ✓ HDF5
- ✓ Plots
- ✓ StatsBase
- ✓ LinearAlgebra, Statistics, Dates, Printf, Random

## Next Steps

1. **For testing only (no data needed)**:
   ```bash
   julia run_tests.jl
   ```

2. **For full benchmarking**:
   - Ensure your HDF5 data file is available (Flt1006_train.h5)
   - Set the `H5_PATH` environment variable
   - Run `julia bench_tl_compare.jl` to generate comparison metrics
   - Results will be saved in `results/` directory

3. **To understand the model**:
   - Read the docstrings in custom_tl.jl for the mathematical model
   - The 9 terms are: 3 permanent + 3 linear-induced + 3 quadratic-induced
   - Each term is normalized for numerical stability

## Troubleshooting

If you encounter issues:
1. Make sure all dependencies are globally installed:
   ```julia
   using Pkg
   Pkg.add(["MagNav", "DataFrames", "CSV", "HDF5", "Plots", "StatsBase"])
   ```

2. If tests fail, check that custom_tl.jl has the Printf import (it's been fixed)

3. For experiment scripts, the artifact system will handle data automatically on first run

## Understanding Julia Artifacts

This project uses **Julia's artifact system** to manage large binary files cleanly:

### How It Works
```julia
# In Julia REPL:
using MagNav
data_path = sgl_2020_train()
println(data_path)  # Shows: ~/.julia/artifacts/<hash>/sgl_2020_train
```

### View Artifacts from Terminal
```bash
# Find the artifact directory
find ~/.julia/artifacts -name "Flt1006_train.h5"

# Or list all artifacts
ls ~/.julia/artifacts
```

### Key Benefits
- ✓ Downloads happen automatically on first use
- ✓ Files are cached for subsequent runs
- ✓ Your repo stays clean (no `git` commits of large binaries)
- ✓ Content-addressed hashing prevents data corruption
- ✓ Works seamlessly across machines
