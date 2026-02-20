#!/usr/bin/env julia
"""
QUICK_START.md

MagNav Benchmark Framework — Complete Quick Start Guide
"""

# 🚀 QUICK START

## 30-Second Setup

```bash
cd /home/lucifer/magnav/magnav-research/magnav-benchmark
julia run_tests.jl
```

✓ Tests pass → You're ready to benchmark!

---

## Run Your First Benchmark

```bash
julia bench_tl_compare.jl
```

That's it. Results are saved in `results/run_<timestamp>/`

---

## What This Project Does

Compares aircraft magnetic compensation algorithms:
- **Custom TL** = Your custom 9-term Tolles-Lawson model
- **Official TL** = MagNav.jl's standard implementation
- **Metrics** = RMS error, DRMS, heading correlation, etc.

---

## Key Features

✓ **No data management** — HDF5 files downloaded automatically  
✓ **Reproducible** — Julia artifact system ensures identical setup everywhere  
✓ **Fast tests** — Unit tests run in <5 seconds (no data needed)  
✓ **Publication-ready** — Multi-segment validation, sensitivity analysis  

---

## All Available Commands

### Testing
```bash
julia run_tests.jl              # Unit tests (10 tests, ~5 seconds)
julia verify_artifact_system.jl # Check data access
```

### Benchmarks
```bash
julia bench_tl_compare.jl       # Main comparison (single segment)
julia bench_multi_segment.jl    # Cross-segment validation
julia bench_sensitivity.jl      # Calibration window length sweep
```

### With Custom Options
```bash
# Custom output directory
RESULTS_DIR=/tmp/results julia bench_tl_compare.jl

# Add experiment notes
RUN_NOTES="testing new model" julia bench_tl_compare.jl
```

---

## Project Structure

```
magnav-benchmark/
├── custom_tl.jl              # 9-term TL model
├── tl_injection.jl           # Compensation injection
├── ekf_wrapper.jl            # EKF benchmarking wrapper
├── nav_metrics.jl            # Performance metrics
├── data_loader.jl            # HDF5 loading utilities
├── experiment_config.jl      # Configuration management
├── segment_utils.jl          # Flight segment utilities
├── bench_plots.jl            # Visualization
│
├── bench_tl_compare.jl       # ← Run this for main benchmark
├── bench_multi_segment.jl    # ← Multi-segment cross-validation
├── bench_sensitivity.jl      # ← Parameter sweep
├── run_tests.jl              # ← Run this first to verify setup
│
├── verify_artifact_system.jl # ← Check data accessibility
├── BenchmarkFramework.jl     # Core module
├── Project.toml              # Dependencies
├── SETUP_SUMMARY.md          # Detailed setup guide
├── ARTIFACT_INTEGRATION.md   # Data path management guide
└── README.md                 # Original framework documentation
```

---

## Data Management (Auto-Handled)

### Available Flights
The artifact system provides 6 flight datasets:
- Flt1002_train.h5 (149 MB)
- Flt1003_train.h5 (120 MB)
- Flt1004_train.h5 (61 MB)
- Flt1005_train.h5 (61 MB)
- **Flt1006_train.h5** (81 MB) ← default
- Flt1007_train.h5 (86 MB)

### Location
Automatically cached at:
```
~/.julia/artifacts/<hash>/sgl_2020_train/Flt1006_train.h5
```

You never need to think about this — it's managed automatically!

### Manual Inspection (if curious)
```bash
# Find exact path
find ~/.julia/artifacts -name "Flt1006_train.h5"

# List all flights
ls ~/.julia/artifacts/*/sgl_2020_train/*.h5

# Check disk usage
du -sh ~/.julia/artifacts
```

---

## Understanding the Code

### Core Classes

**CustomTLModel**
```julia
struct CustomTLModel
    coef::Vector{Float64}        # 9 fitted coefficients
    residual_rms::Float64        # RMS of residuals (nT)
    var_reduction::Float64       # Variance reduction (0-1)
    heading_corr_post::Float64   # Heading corr. after compensation
end
```

**BenchmarkData**
```julia
struct BenchmarkData
    xyz::XYZ20         # Flight trajectory & sensor data
    mapS::MapS         # Magnetic anomaly map
    flight_id::Symbol  # e.g., :Flt1006
    h5_path::String    # Source HDF5 file
end
```

### Key Functions

| Function | Purpose |
|----------|---------|
| `build_A9()` | Build 9-term TL design matrix |
| `fit_custom_tl!()` | Fit custom TL model (OLS/ridge) |
| `compensate_flight_data()` | Apply compensation to raw data |
| `compute_nav_metrics()` | Calculate navigation errors |
| `load_benchmark_data()` | Load HDF5 → structured data |

---

## Benchmark Output

Results saved to `results/run_<YYYYMMDD_HHMMSS>/`:

```
results/
└── run_20260217_150623/
    ├── config.txt              # Experiment configuration
    ├── figures/
    │   ├── trajectory_comparison.png
    │   ├── heading_correlation.png
    │   ├── error_distribution.png
    │   └── residual_timeseries.png
    └── tables/
        ├── metrics_summary.csv
        └── segment_comparison.csv
```

---

## Troubleshooting

### Tests fail?
```bash
# Check Julia is 1.10+
julia --version

# Verify dependencies
julia -e 'using MagNav; using DataFrames; using StatsBase; println("OK")'
```

### Can't find artifact data?
```bash
# Verify artifact system is working
julia verify_artifact_system.jl
```

### Benchmark script fails?
Check the error message — usually one of:
1. Missing dependency → `using Pkg; Pkg.add("DependencyName")`
2. Data not found → Run `verify_artifact_system.jl`
3. Insufficient disk space → `du -sh ~/.julia/artifacts`

---

## Next Steps

1. **Run tests** (verify setup works):
   ```bash
   julia run_tests.jl
   ```

2. **Run main benchmark** (compare methods):
   ```bash
   julia bench_tl_compare.jl
   ```

3. **Explore results**:
   - Check generated plots in `results/run_*/figures/`
   - Check metric tables in `results/run_*/tables/`

4. **Modify and experiment**:
   - Edit `custom_tl.jl` to try different TL models
   - Change `experiment_config.jl` for different parameters
   - Run `bench_multi_segment.jl` for cross-segment validation

---

## Learning Resources

### Inside This Project
- [SETUP_SUMMARY.md](SETUP_SUMMARY.md) — Detailed setup guide
- [ARTIFACT_INTEGRATION.md](ARTIFACT_INTEGRATION.md) — Data management internals
- [README.md](README.md) — Original framework documentation
- [custom_tl.jl](custom_tl.jl) — Mathematical model with docstrings

### External
- **MagNav.jl** — Check MagNav documentation for official TL details
- **Julia Artifacts** — https://pkgdocs.julialang.org/v1/artifacts/

---

## Tips & Tricks

### Faster Iterations During Development
```julia
# In Julia REPL from project directory:
include("custom_tl.jl")
data = synthetic_tl_data(n=5000)
model = CustomTLModel()
fit_custom_tl!(model, data.Bx, data.By, data.Bz, data.mag, data.heading; verbose=true)
```

### Compare Different Flight Data
Edit benchmark script to change:
```julia
const H5_PATH = joinpath(H5_DIR, "Flt1003_train.h5")  # Try different flight
```

### Profile a Benchmark Run
```bash
julia --profile=cpu bench_tl_compare.jl
```

### Parallel Benchmarking (if you have multiple cores)
```bash
julia -p auto bench_multi_segment.jl  # Use all CPU cores
```

---

## Summary

**This benchmark framework is now fully automated:**

```bash
# Tests (no data)
julia run_tests.jl           # ✓ Works immediately

# Benchmarks (data auto-downloaded)
julia bench_tl_compare.jl    # ✓ Works immediately
```

No manual configuration needed. All HDF5 data is managed by Julia's artifact system.

**Questions?** Check [ARTIFACT_INTEGRATION.md](ARTIFACT_INTEGRATION.md) or run:
```bash
julia verify_artifact_system.jl
```

Happy benchmarking! 🚀
