# 📌 REFERENCE CARD

## One-Liner Commands

```bash
# Run tests (no data needed, ~5 seconds)
julia run_tests.jl

# Run main benchmark (auto-downloads data)
julia bench_tl_compare.jl

# Run all experiments
julia bench_tl_compare.jl && julia bench_multi_segment.jl && julia bench_sensitivity.jl

# Verify data access
julia verify_artifact_system.jl
```

## Where is my data?

```bash
~/.julia/artifacts/<hash>/sgl_2020_train/Flt1006_train.h5
```

**Don't worry about it — Julia manages it automatically!**

## Documentation Map

| I want to... | Read this |
|---|---|
| Get started in 30 seconds | [QUICK_START.md](QUICK_START.md) |
| Understand artifact system | [ARTIFACT_INTEGRATION.md](ARTIFACT_INTEGRATION.md) |
| Full setup guide | [SETUP_SUMMARY.md](SETUP_SUMMARY.md) |
| See what changed | This file (INTEGRATION_COMPLETE.md) |
| Learn the math | [custom_tl.jl](custom_tl.jl) |
| Understand benchmarks | [README.md](README.md) |

## Environment Variables (Optional)

```bash
# Custom output directory
RESULTS_DIR=/tmp/results julia bench_tl_compare.jl

# Add notes to results
RUN_NOTES="my test run" julia bench_tl_compare.jl

# Both together
RESULTS_DIR=/tmp/results RUN_NOTES="test" julia bench_tl_compare.jl
```

## Installed Test Status

```
✓ Julia 1.10.0
✓ MagNav installed
✓ DataFrames installed
✓ StatsBase installed
✓ All dependencies ready
✓ Tests: 10/10 passing
✓ Artifact system: Working
✓ Flight data: All 6 available
```

## If You Need Help

1. **Quick overview:** `cat QUICK_START.md`
2. **Something broken:** `julia verify_artifact_system.jl`
3. **Want to understand:** `cat ARTIFACT_INTEGRATION.md`
4. **Full details:** `cat SETUP_SUMMARY.md`

## Timeline

| When | Action | Result |
|------|--------|--------|
| Now | `julia run_tests.jl` | ✓ Tests pass |
| Next | `julia bench_tl_compare.jl` | ✓ Benchmark runs |
| Minutes | Check results | ✓ Plots generated |
| Later | Modify code | ✓ Re-run benchmarks |

## Key Files in Project

```
Core Algorithms
├── custom_tl.jl ..................... 9-term TL model
├── ekf_wrapper.jl ................... EKF navigation filter
└── nav_metrics.jl ................... Error computations

Utilities
├── data_loader.jl ................... Load HDF5 → data structures
├── segment_utils.jl ................. Flight segment handling
├── experiment_config.jl ............. Configuration management
└── tl_injection.jl .................. Inject compensation

Visualization
└── bench_plots.jl ................... Generate figures

Benchmarks (RUN THESE)
├── bench_tl_compare.jl .............. Main benchmark
├── bench_multi_segment.jl ........... Multi-segment validation
├── bench_sensitivity.jl ............. Parameter sweep
├── run_tests.jl ..................... Unit tests
└── verify_artifact_system.jl ........ Data verification

Documentation
├── QUICK_START.md ................... Start here
├── ARTIFACT_INTEGRATION.md .......... Data management details
├── SETUP_SUMMARY.md ................. Full setup reference
├── INTEGRATION_COMPLETE.md .......... This summary
└── README.md ........................ Original framework docs
```

## Results Directory Structure

After running `julia bench_tl_compare.jl`:

```
results/
└── run_20260217_150623/
    ├── config.txt
    ├── figures/
    │   ├── trajectory_comparison.png
    │   ├── heading_correlation.png
    │   └── ...
    └── tables/
        ├── metrics_summary.csv
        └── ...
```

## Quick Debugging

```bash
# Check if Julia works
julia --version

# Check if MagNav works
julia -e "using MagNav; println(\"OK\")"

# Check if artifacts work
julia verify_artifact_system.jl

# Check if tests pass
julia run_tests.jl

# Check current directory
pwd

# Check disk space
du -sh ~/.julia/artifacts
```

## Pro Tips

1. **First time running?** → Expect ~1 min for data download
2. **Subsequent runs?** → Much faster (cached data)
3. **Different flight?** → Edit the H5_PATH line in benchmark script
4. **Parallel processing?** → Use `julia -p auto bench_*.jl`
5. **Working directory?** → Should be the magnav-benchmark folder

## What Happens When You Run `julia bench_tl_compare.jl`

```
1. Script loads → checks dependencies
2. Calls sgl_2020_train() → downloads artifact (if needed)
3. Opens Flt1006_train.h5 → loads flight data
4. Splits data → calibration + evaluation segments
5. Fits custom TL → optimizes 9 coefficients
6. Compensates magnetic data → applies correction
7. Runs EKF → navigation simulation
8. Computes metrics → RMS, DRMS, etc.
9. Saves results → plots + CSV files
10. Done! → Results in results/run_<timestamp>/
```

## FAQ

**Q: Why does it take time on first run?**  
A: Downloading and caching the 600MB+ artifact directory.

**Q: Can I use my own data?**  
A: Yes, modify the `load_benchmark_data()` function in data_loader.jl.

**Q: How do I use different flights?**  
A: Change `H5_PATH = joinpath(H5_DIR, "Flt1003_train.h5")` in benchmark script.

**Q: Is data stored in my repo?**  
A: No! It's in ~/.julia/artifacts/. Your repo stays clean.

**Q: What if artifact system fails?**  
A: Run `julia verify_artifact_system.jl` to diagnose.

---

**Status:** ✓ Ready to Go  
**Last Check:** Tests passing (10/10)  
**Next Step:** `julia bench_tl_compare.jl`
