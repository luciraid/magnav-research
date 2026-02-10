# Diurnal Magnetic Field Correction – Flt1002

## Action
- Loaded raw scalar magnetometer data (`mag_1_c`)
- Loaded diurnal reference signal
- Applied diurnal correction by subtraction

## Result
- Diurnal-corrected signal shows reduced low-frequency temporal variation
- Remaining structure is dominated by aircraft interference and geology

## Output
- figures/Flt1002_raw_vs_diurnal_corrected.png

## Relevance to MagNav
Diurnal correction is a prerequisite step before aircraft magnetic
interference compensation (Tolles–Lawson) and magnetic anomaly matching.
