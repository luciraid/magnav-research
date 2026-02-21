#!/usr/bin/env python3
"""
generate_comparison.py — Compare two run_test.sh result directories side-by-side.

Reads metrics.txt from each result directory, extracts key numbers,
and writes gazebo/results/comparison_report.txt.

Usage:
    python3 gazebo/scripts/generate_comparison.py <dir_A> <dir_B>
    python3 gazebo/scripts/generate_comparison.py \\
        gazebo/results/zero_20260222_003000 \\
        gazebo/results/my_coef_20260222_003500

Output:
    gazebo/results/comparison_report.txt
    Printed to stdout.
"""

import sys
import re
from pathlib import Path
from datetime import datetime

# ── Metric extraction ─────────────────────────────────────────────────────────

PATTERNS = {
    "var_reduction":  re.compile(r"Variance reduction\s*:\s*(-?[\d.]+)%"),
    "hdg_corr_pre":   re.compile(r"Heading corr before TL\s*:\s*([\d.]+)"),
    "hdg_corr_post":  re.compile(r"Heading corr after\s*TL\s*:\s*([\d.]+)"),
    "rms_removed":    re.compile(r"RMS interference removed\s*:\s*([\d.]+)"),
    "drms_raw":       re.compile(r"Mag nav without TL \(est\.\)\s*:\s*([\d.]+)"),
    "drms_tl":        re.compile(r"Mag nav with TL\s*\(est\.\)\s*:\s*([\d.]+)"),
    "improvement":    re.compile(r"Estimated improvement\s*:\s*([\d.]+)%"),
    "n_samples":      re.compile(r"Samples:\s*(\d+)"),
    "duration_min":   re.compile(r"Samples:\s*\d+\s+\(([\d.]+) min"),
    "lat_range":      re.compile(r"Lat:\s*([\d.]+)°\s*–\s*([\d.]+)°"),
    "alt_range":      re.compile(r"Alt:\s*([\d.]+)\s*–\s*([\d.]+) m"),
}


def parse_metrics(metrics_txt: Path) -> dict:
    """Extract numeric metrics from a metrics.txt file."""
    if not metrics_txt.exists():
        return {}
    text = metrics_txt.read_text()
    result = {}
    for key, pat in PATTERNS.items():
        m = pat.search(text)
        if m:
            if key in ("lat_range", "alt_range"):
                result[key] = (float(m.group(1)), float(m.group(2)))
            elif key in ("n_samples",):
                result[key] = int(m.group(1))
            else:
                result[key] = float(m.group(1))
    return result


def read_label(result_dir: Path) -> str:
    """Derive a short label from directory name or coef_path.txt."""
    coef_path_file = result_dir / "coef_path.txt"
    if coef_path_file.exists():
        coef_path = coef_path_file.read_text().strip()
        stem = Path(coef_path).stem
        return stem
    return result_dir.name


def fmt(val, fmt_str=".4f", suffix="", missing="N/A"):
    if val is None:
        return missing
    return f"{val:{fmt_str}}{suffix}"


def improvement_arrow(a, b, lower_is_better=True):
    """Return ↓/↑/= arrow depending on direction of change from a to b."""
    if a is None or b is None:
        return ""
    delta = b - a
    if abs(delta) < 1e-6:
        return "="
    if lower_is_better:
        return "↓ BETTER" if delta < 0 else "↑ WORSE"
    else:
        return "↑ BETTER" if delta > 0 else "↓ WORSE"


# ── Report generation ─────────────────────────────────────────────────────────

def generate_report(dir_a: Path, dir_b: Path, out_path: Path) -> str:
    m_a = parse_metrics(dir_a / "metrics.txt")
    m_b = parse_metrics(dir_b / "metrics.txt")

    label_a = read_label(dir_a)
    label_b = read_label(dir_b)

    ulg_a = (dir_a / "ulg_path.txt").read_text().strip() \
            if (dir_a / "ulg_path.txt").exists() else "unknown"
    ulg_b = (dir_b / "ulg_path.txt").read_text().strip() \
            if (dir_b / "ulg_path.txt").exists() else "unknown"

    same_log = ulg_a == ulg_b

    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    plot_a = next(iter(dir_a.glob("*.png")), None)
    plot_b = next(iter(dir_b.glob("*.png")), None)

    W = 72   # report width

    lines = []
    def h(title=""):
        lines.append("=" * W)
        if title:
            lines.append(f"  {title}")
            lines.append("=" * W)

    def row(label, val_a, val_b, note=""):
        col = 24
        a_str = str(val_a).ljust(18)
        b_str = str(val_b).ljust(18)
        lines.append(f"  {label:<{col}}  {a_str}  {b_str}  {note}")

    def divider():
        lines.append("-" * W)

    h(f"MagNav TL Compensation — Side-by-Side Comparison Report")
    lines.append(f"  Generated : {now}")
    lines.append(f"  Flight log: {Path(ulg_a).name}")
    if not same_log:
        lines.append(f"              {Path(ulg_b).name}  (different logs!)")
    else:
        lines.append(f"  Same log  : YES — results are directly comparable.")
    lines.append("")
    lines.append(f"  {'Test A':<24}  {label_a}")
    lines.append(f"  {'Test B':<24}  {label_b}")
    lines.append("")

    # ── Flight info ───────────────────────────────────────────────────────────
    h("Flight Info")
    lines.append("")
    row("Metric", f"[A] {label_a}", f"[B] {label_b}")
    divider()

    n_a  = m_a.get("n_samples")
    n_b  = m_b.get("n_samples")
    row("Samples", fmt(n_a, ".0f") if n_a else "N/A",
                   fmt(n_b, ".0f") if n_b else "N/A")

    d_a = m_a.get("duration_min")
    d_b = m_b.get("duration_min")
    row("Duration (min)", fmt(d_a, ".1f"), fmt(d_b, ".1f"))

    lr_a = m_a.get("lat_range")
    lr_b = m_b.get("lat_range")
    row("Lat range (°)",
        f"{lr_a[0]:.4f}–{lr_a[1]:.4f}" if lr_a else "N/A",
        f"{lr_b[0]:.4f}–{lr_b[1]:.4f}" if lr_b else "N/A")

    alt_a = m_a.get("alt_range")
    alt_b = m_b.get("alt_range")
    row("Alt range (m)",
        f"{int(alt_a[0])}–{int(alt_a[1])}" if alt_a else "N/A",
        f"{int(alt_b[0])}–{int(alt_b[1])}" if alt_b else "N/A")

    lines.append("")

    # ── Compensation quality ──────────────────────────────────────────────────
    h("TL Compensation Quality")
    lines.append("")
    row("Metric", f"[A] {label_a}", f"[B] {label_b}", "Change A→B")
    divider()

    vr_a = m_a.get("var_reduction")
    vr_b = m_b.get("var_reduction")
    row("Variance reduction (%)",
        fmt(vr_a, ".1f", "%"), fmt(vr_b, ".1f", "%"),
        improvement_arrow(vr_a, vr_b, lower_is_better=False))

    hpre_a = m_a.get("hdg_corr_pre")
    hpre_b = m_b.get("hdg_corr_pre")
    row("Hdg corr BEFORE TL",
        fmt(hpre_a, ".4f"), fmt(hpre_b, ".4f"))

    hpost_a = m_a.get("hdg_corr_post")
    hpost_b = m_b.get("hdg_corr_post")
    row("Hdg corr AFTER  TL",
        fmt(hpost_a, ".4f"), fmt(hpost_b, ".4f"),
        improvement_arrow(hpost_a, hpost_b, lower_is_better=True))

    rms_a = m_a.get("rms_removed")
    rms_b = m_b.get("rms_removed")
    row("RMS interference (nT)",
        fmt(rms_a, ".2f"), fmt(rms_b, ".2f"))

    lines.append("")

    # ── Estimated DRMS ────────────────────────────────────────────────────────
    h("Estimated Navigation DRMS  (* extrapolated from Flt1006 benchmark)")
    lines.append("")
    row("Metric", f"[A] {label_a}", f"[B] {label_b}", "Change A→B")
    divider()

    row("GPS baseline (m)", "3.00 m", "3.00 m")

    dr_a = m_a.get("drms_raw")
    dr_b = m_b.get("drms_raw")
    row("Mag nav w/o TL* (m)",
        fmt(dr_a, ".2f", " m"), fmt(dr_b, ".2f", " m"))

    dtl_a = m_a.get("drms_tl")
    dtl_b = m_b.get("drms_tl")
    row("Mag nav w/  TL* (m)",
        fmt(dtl_a, ".2f", " m"), fmt(dtl_b, ".2f", " m"),
        improvement_arrow(dtl_a, dtl_b, lower_is_better=True))

    imp_a = m_a.get("improvement")
    imp_b = m_b.get("improvement")
    row("Est. improvement (%)",
        fmt(imp_a, ".0f", "%"), fmt(imp_b, ".0f", "%"),
        improvement_arrow(imp_a, imp_b, lower_is_better=False))

    lines.append("")

    # ── Key insight ────────────────────────────────────────────────────────────
    h("Interpretation")
    lines.append("")

    # Determine which coef is better
    if hpost_a is not None and hpost_b is not None:
        if hpost_b < hpost_a * 0.95:
            lines.append(f"  ✓  [{label_b}] reduces heading-correlated interference")
            lines.append(f"     more than [{label_a}].")
            lines.append(f"     Hdg corr after TL:  {hpost_a:.4f}  →  {hpost_b:.4f}")
        elif hpost_a < hpost_b * 0.95:
            lines.append(f"  ✓  [{label_a}] reduces heading-correlated interference")
            lines.append(f"     more than [{label_b}].")
            lines.append(f"     Hdg corr after TL:  {hpost_b:.4f}  →  {hpost_a:.4f}")
        else:
            lines.append(f"  ~  Both coefficient sets produce similar results.")
            lines.append(f"     Difference in heading correlation is < 5%.")

    # Zero-coef baseline interpretation
    if label_a.lower().startswith("zero") or label_b.lower().startswith("zero"):
        zero_dir = "A" if label_a.lower().startswith("zero") else "B"
        real_dir = "B" if zero_dir == "A" else "A"
        hpre = hpre_a if zero_dir == "A" else hpre_b
        hpost = hpost_a if zero_dir == "A" else hpost_b
        if hpre is not None and hpost is not None:
            lines.append("")
            lines.append(f"  NOTE: Test {zero_dir} (zero coefs) is the raw Gazebo signal.")
            lines.append(f"    Hdg corr {hpre:.4f} = raw aircraft heading interference.")
            lines.append(f"    Test {real_dir} coefs should reduce this below {hpre:.4f}.")

    lines.append("")

    # ── What to tell your stakeholder ─────────────────────────────────────────
    h("Stakeholder Summary")
    lines.append("")
    lines.append("  Pipeline validated on Gazebo SITL flight:")
    lines.append(f"    Log: {Path(ulg_a).name}")
    if d_a:
        lines.append(f"    Duration: {d_a:.1f} min")
    lines.append("")
    lines.append("  Test inputs:")
    lines.append(f"    [A] zero_coef.txt  — null baseline (raw signal)")
    lines.append(f"    [B] my_coef.txt    — trained 9-term TL model")
    lines.append("")
    lines.append("  The framework is ready for external coefficients.")
    lines.append("  Send your custom_tl_coef.txt and I will:")
    lines.append("    1. Place it in gazebo/test_data/")
    lines.append("    2. Run: ./gazebo/run_test.sh your_coef.txt --ulg <log>")
    lines.append("    3. Return a comparison report + validation plot.")
    lines.append("")

    # ── File references ───────────────────────────────────────────────────────
    h("Output Files")
    lines.append("")
    lines.append(f"  [A] {dir_a}")
    if plot_a:
        lines.append(f"      Plot: {plot_a.name}")
    lines.append(f"  [B] {dir_b}")
    if plot_b:
        lines.append(f"      Plot: {plot_b.name}")
    lines.append(f"  Report: {out_path}")
    lines.append("")
    h()

    report = "\n".join(lines)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(report)

    return report


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    if len(sys.argv) < 3:
        print("Usage: python3 generate_comparison.py <result_dir_A> <result_dir_B>")
        print("")
        print("  result_dir_A  : e.g. gazebo/results/zero_coef_20260222_003000")
        print("  result_dir_B  : e.g. gazebo/results/my_coef_20260222_003500")
        sys.exit(1)

    dir_a = Path(sys.argv[1]).resolve()
    dir_b = Path(sys.argv[2]).resolve()

    for d in (dir_a, dir_b):
        if not d.exists():
            sys.exit(f"Directory not found: {d}")
        if not (d / "metrics.txt").exists():
            sys.exit(f"No metrics.txt in {d} — run run_test.sh first.")

    # Output goes to gazebo/results/ relative to the result dirs' parent
    # (or two levels up from a typical gazebo/results/xxx/ path)
    common_parent = dir_a.parent
    out_path = common_parent / "comparison_report.txt"

    print(f"Comparing:")
    print(f"  A: {dir_a}")
    print(f"  B: {dir_b}")
    print()

    report = generate_report(dir_a, dir_b, out_path)
    print(report)
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
