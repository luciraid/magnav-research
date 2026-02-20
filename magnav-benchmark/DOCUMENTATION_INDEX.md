# 📚 Documentation Index

## 🚀 START HERE

### For the Impatient
**File:** [START_HERE.md](START_HERE.md)  
**Time:** 2 minutes  
**Contains:** Project status, what to run next

### For the Curious  
**File:** [QUICK_START.md](QUICK_START.md)  
**Time:** 5 minutes  
**Contains:** 30-second setup, key concepts, all commands

### For the Reference Lookups
**File:** [REFERENCE_CARD.md](REFERENCE_CARD.md)  
**Time:** 2 minutes  
**Contains:** One-liners, FAQ, quick debugging

---

## 📖 DETAILED GUIDES

### Understanding the Artifact System
**File:** [ARTIFACT_INTEGRATION.md](ARTIFACT_INTEGRATION.md)  
**Time:** 15 minutes  
**Contains:** How data is managed, benefits, internals

### Complete Setup Reference
**File:** [SETUP_SUMMARY.md](SETUP_SUMMARY.md)  
**Time:** 10 minutes  
**Contains:** Classes, functions, dependencies, troubleshooting

### What Changed in This Session
**File:** [INTEGRATION_COMPLETE.md](INTEGRATION_COMPLETE.md)  
**Time:** 10 minutes  
**Contains:** Before/after, changes made, verification checklist

---

## 🔧 CODE DOCUMENTATION

### Mathematical Model
**File:** [custom_tl.jl](custom_tl.jl)  
**Focus:** 9-term Tolles-Lawson model, build_A9(), fit_custom_tl!()  
**Read:** Always check docstrings!

### Original Framework Docs
**File:** [README.md](README.md)  
**Contains:** Pipeline overview, directory structure, experimental design

---

## 🎯 QUICK NAVIGATION

### I Want To...

| Goal | File | Time |
|------|------|------|
| **Get running ASAP** | [START_HERE.md](START_HERE.md) | 2 min |
| **Understand what I have** | [QUICK_START.md](QUICK_START.md) | 5 min |
| **Look up a command** | [REFERENCE_CARD.md](REFERENCE_CARD.md) | 2 min |
| **Deep dive on data** | [ARTIFACT_INTEGRATION.md](ARTIFACT_INTEGRATION.md) | 15 min |
| **Full reference** | [SETUP_SUMMARY.md](SETUP_SUMMARY.md) | 10 min |
| **See what changed** | [INTEGRATION_COMPLETE.md](INTEGRATION_COMPLETE.md) | 10 min |
| **Learn the math** | [custom_tl.jl](custom_tl.jl) docstrings | 20 min |

---

## 🗂️ FILE ORGANIZATION

### Documentation Files (This Session)
```
START_HERE.md           ← You are here (entry point)
QUICK_START.md          ← 30-second getting started
REFERENCE_CARD.md       ← Command cheat sheet
SETUP_SUMMARY.md        ← Complete reference guide
ARTIFACT_INTEGRATION.md ← Data management deep dive
INTEGRATION_COMPLETE.md ← Session summary
DOCUMENTATION_INDEX.md  ← This file
```

### Core Files (Always Check Docstrings)
```
custom_tl.jl            ← 9-term TL model
ekf_wrapper.jl          ← Navigation filter
nav_metrics.jl          ← Error computation
data_loader.jl          ← HDF5 loading
segment_utils.jl        ← Flight segments
experiment_config.jl    ← Configuration
```

### Benchmark Scripts
```
bench_tl_compare.jl     ← Main benchmark (run this)
bench_multi_segment.jl  ← Multi-segment validation
bench_sensitivity.jl    ← Parameter sweep
run_tests.jl            ← Unit tests
verify_artifact_system.jl ← Data verification
```

---

## 📋 DOCUMENTATION MATRIX

```
                    | TIME | TYPE       | DETAIL LEVEL
            --------+------+------------|---------------
START_HERE  | 2m   | Guide      | Overview
QUICK_START | 5m   | Guide      | Practical
REFERENCE   | 2m   | Reference  | Commands
SETUP_SUMMARY | 10m | Reference | Complete
ARTIFACT    | 15m  | Deep Dive  | Technical
INTEGRATION | 10m  | Summary    | Session
```

---

## 🎓 RECOMMENDED READING ORDER

### For Complete Beginners
1. [START_HERE.md](START_HERE.md) — Understand what you have
2. [QUICK_START.md](QUICK_START.md) — Learn how to use it
3. [REFERENCE_CARD.md](REFERENCE_CARD.md) — Get commands reference
4. Run `julia run_tests.jl` — Verify it works
5. Run `julia bench_tl_compare.jl` — See it in action
6. [custom_tl.jl](custom_tl.jl) — Learn the algorithms

### For Developers
1. [INTEGRATION_COMPLETE.md](INTEGRATION_COMPLETE.md) — See what changed
2. [ARTIFACT_INTEGRATION.md](ARTIFACT_INTEGRATION.md) — Understand architecture
3. [SETUP_SUMMARY.md](SETUP_SUMMARY.md) — Full reference
4. [custom_tl.jl](custom_tl.jl) — Study the code
5. Modify and experiment

### For Quick Reference
1. [REFERENCE_CARD.md](REFERENCE_CARD.md) — Get the command
2. Run it
3. Check results in `results/`

---

## 💡 WHICH FILE TO READ?

### Question: "How do I get started?"
→ [START_HERE.md](START_HERE.md) or [QUICK_START.md](QUICK_START.md)

### Question: "What command do I need?"
→ [REFERENCE_CARD.md](REFERENCE_CARD.md)

### Question: "How does the artifact system work?"
→ [ARTIFACT_INTEGRATION.md](ARTIFACT_INTEGRATION.md)

### Question: "What was changed?"
→ [INTEGRATION_COMPLETE.md](INTEGRATION_COMPLETE.md)

### Question: "Complete setup guide?"
→ [SETUP_SUMMARY.md](SETUP_SUMMARY.md)

### Question: "How does the math work?"
→ [custom_tl.jl](custom_tl.jl) and [README.md](README.md)

---

## ✅ VERIFICATION

### You Know Everything You Need If You Can Answer:

- [ ] What command runs the tests?
- [ ] What command runs the main benchmark?
- [ ] Where is the HDF5 data stored?
- [ ] What does the artifact system do?
- [ ] How do you provide a custom results directory?

### Answers:

1. `julia run_tests.jl`
2. `julia bench_tl_compare.jl`
3. `~/.julia/artifacts/<hash>/sgl_2020_train/`
4. Automatically downloads and caches data
5. `RESULTS_DIR=/path julia bench_tl_compare.jl`

→ If you can answer these, you're ready!

---

## 🚀 QUICK START FOR DIFFERENT PERSONAS

### Impatient Engineer
```bash
julia run_tests.jl
julia bench_tl_compare.jl
# Check results in results/
```

### Careful Researcher
```bash
# 1. Read the overview
cat START_HERE.md

# 2. Verify everything works
julia run_tests.jl
julia verify_artifact_system.jl

# 3. Run benchmarks
julia bench_tl_compare.jl
julia bench_multi_segment.jl

# 4. Analyze results
ls -R results/
```

### Academic/Publisher
```bash
# 1. Understand methodology
cat README.md
cat ARTIFACT_INTEGRATION.md

# 2. Reproduce results
julia bench_tl_compare.jl
julia bench_multi_segment.jl
julia bench_sensitivity.jl

# 3. Generate figures
# Figures auto-generated in results/

# 4. Extract metrics
cat results/run_*/tables/metrics_summary.csv
```

---

## 📞 TROUBLESHOOTING QUICK LINKS

| Problem | Solution |
|---------|----------|
| "How do I start?" | Read [START_HERE.md](START_HERE.md) |
| "I need a command" | Check [REFERENCE_CARD.md](REFERENCE_CARD.md) |
| "Data not found" | Run `julia verify_artifact_system.jl` |
| "Need full setup guide" | Read [SETUP_SUMMARY.md](SETUP_SUMMARY.md) |
| "Want technical details" | Read [ARTIFACT_INTEGRATION.md](ARTIFACT_INTEGRATION.md) |
| "Tests failing" | Check [INTEGRATION_COMPLETE.md](INTEGRATION_COMPLETE.md) verification |

---

## 📊 DOCUMENTATION STATISTICS

```
Total Documentation:  ~10,000 words
Total Code Comments:  ~2,000 lines
Guides Provided:      7 files
Immediate Start Time: 2 minutes
Full Understanding:   30 minutes
Ready to Benchmark:   5 minutes
```

---

## 🎯 READING TIME ESTIMATES

| Document | Time | Best For |
|---|---:|---|
| START_HERE | 2 min | Getting momentum |
| QUICK_START | 5 min | Understanding project |
| REFERENCE_CARD | 2 min | Lookup help |
| SETUP_SUMMARY | 10 min | Complete reference |
| ARTIFACT_INTEGRATION | 15 min | Technical depth |
| INTEGRATION_COMPLETE | 10 min | Session summary |
| **TOTAL** | **44 min** | **Full mastery** |

Or just run `julia bench_tl_compare.jl` and learn by doing! 🚀

---

## 📌 KEY TAKEAWAYS

1. **Start immediately:** [START_HERE.md](START_HERE.md)
2. **One command:** `julia bench_tl_compare.jl`
3. **Data auto-managed:** No manual setup needed
4. **Tests verify setup:** `julia run_tests.jl`
5. **Full docs available:** Check this index

---

**Last Updated:** February 17, 2026  
**Files:** 7 documentation + code with docstrings  
**Ready to Use:** Yes ✓
