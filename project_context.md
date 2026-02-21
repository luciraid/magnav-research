# MagNav Research – Project Context

## Environment
- OS: Ubuntu 24.04
- Julia: 1.10
- IDE: VS Code
- MagNav cloned locally:
  ~/magnav/MagNav.jl
- Custom research repo:
  ~/magnav/magnav-research

---

## Dataset
- Source: DAF-MIT AIA Open Flight Data (SGL 2020)
- Accessed via: sgl_2020_train()
- Flight currently used: Flt1006_train.h5
- Map used: ottawa_area_maps artifact

---

## Objective

Build a reproducible MagNav benchmarking framework comparing:

1. Official MagNav Tolles–Lawson compensation
2. Custom 9-term TL implementation

Evaluation level: Navigation performance (not just signal quality).

---

## Current Architecture

### Data Loading
- sgl_2020_train()
- h5_to_df(...)
- get_XYZ(:Flt1006, df; silent=true)

### Official TL
- create_TL_coef(...)
- create_TL_A(...)
- detrend(A * TL_coef; mean_only=true)

### Navigation
- run_filt(...)
- EKF
- Map interpolation via get_map_val(...)

---

## Current Performance (Official TL)

- INS DRMS ≈ 137 m
- MagNav DRMS ≈ 31 m
- CRLB ≈ 13 m

Navigation error remains bounded (~10–20 m range)
INS-only drifts to ~200 m in ~14 min

---

## Custom TL Status

- Implemented standalone 9-term TL
- Validated:
  - ~18% variance reduction
  - Heading correlation reduced ~0.08 → ~0.001
- Eddy terms unstable without regularization
- Currently not yet injected into MagNav pipeline

---

## Benchmark Script

Location:
~/magnav/magnav-research/src/benchmark_tl.jl

Capabilities:
- Runs official TL
- Runs EKF
- Computes:
  - RMS North
  - RMS East
  - DRMS
  - Final position error
- Placeholder for custom TL insertion

---

## Immediate Next Goal

Inject custom TL compensation into benchmark script and compare:

- DRMS (official vs custom)
- RMS North/East
- Final position error
- Heading correlation
- Error vs time plots

No Pluto.
No modifying MagNav internal package code.
Fully reproducible from project root.

---

## Long-Term Goals

- Determine if custom TL improves navigation performance
- Quantify sensitivity to calibration window
- Analyze observability vs flight segment
- Evaluate robustness to INS noise
- Prepare technical summary for stakeholder

---

## Constraints

- Clean modular architecture
- No notebook-only context
- No hidden state
- Fully script-based execution
- Results reproducible

---

## Open Questions

- Does custom TL improve DRMS?
- Is navigation improvement statistically significant?
- Does compensation distort anomaly band?
- How sensitive is performance to calibration window length?

---

End of context file.
